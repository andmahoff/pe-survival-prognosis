"""p2_nn.py
Small neural networks in numpy, imported by 10_deep.py (per-modality neural
learners) and 18j_arch_neural.py (intermediate fusion).

The network takes one or more input blocks, one per modality. Each block has
its own encoder (dense layer + ReLU); the encodings are concatenated and
passed through dropout, a trunk (dense layer + ReLU + dropout) and a linear
head. A single block is a plain two-layer MLP. Multimodal dropout (training
only) drops each block's encoding with probability md, always keeping at
least one block, and rescales the kept encodings by 1 / (1 - md).

Each loss returns the loss and its gradient with respect to the head:
  cox       negative Cox partial log-likelihood (Breslow ties), the
            DeepSurv loss (Katzman et al. 2018)
  bce       binary cross-entropy for several outputs at once (the MLP
            binary twin: one output per horizon)
  hazard    discrete-time logistic hazard (Kvamme and Borgan 2021):
            one output per time interval, binary cross-entropy over
            the intervals each patient was at risk in
  deephit   DeepHit for competing risks (Lee et al. 2018): a softmax
            over every (cause, interval) cell plus one 'no event'
            cell; loss = alpha * negative log-likelihood + (1 -
            alpha) * the cause-specific ranking loss with scale
            sigma (alpha and sigma can travel in the target
            dictionary)

Training uses Adam, mini-batches and L2 weight decay, with early stopping on
a validation set (best weights kept) or a fixed number of epochs.

Checks:  grad_check() compares every analytic gradient with finite
         differences.
Usage:   import p2_nn as nn   (as in 10_deep.py and 18j_arch_neural.py)
"""
import copy
import numpy as np


def relu(z):
    return np.maximum(z, 0.0)


class Net:
    def __init__(self, blocks, enc, trunk, out, dropout=0.1, md=0.0,
                 rng=None):
        self.blocks = [np.asarray(b, int) for b in blocks]
        self.enc = list(enc)
        self.dropout = float(dropout)
        self.md = float(md)
        rng = rng if rng is not None else np.random.RandomState(0)
        P = {}
        for i, (b, e) in enumerate(zip(self.blocks, self.enc)):
            P["W%d" % i] = rng.randn(len(b), e) * np.sqrt(2.0 / len(b))
            P["c%d" % i] = np.zeros(e)
        h0 = int(sum(self.enc))
        P["V"] = rng.randn(h0, trunk) * np.sqrt(2.0 / h0)
        P["v"] = np.zeros(trunk)
        P["U"] = rng.randn(trunk, out) * np.sqrt(1.0 / trunk)
        P["u"] = np.zeros(out)
        self.P = P

    def forward(self, X, train=False, rng=None):
        P = self.P
        cache = {"X": X, "Z": [], "A": [], "M": []}
        n = X.shape[0]
        nb = len(self.blocks)
        keep = np.ones((n, nb))
        if train and self.md > 0 and nb > 1:
            keep = (rng.rand(n, nb) >= self.md).astype(float)
            none = keep.sum(1) == 0
            if none.any():
                keep[none, rng.randint(nb, size=none.sum())] = 1.0
            keep = keep / (1.0 - self.md)
        enc_out = []
        for i, b in enumerate(self.blocks):
            Z = X[:, b] @ P["W%d" % i] + P["c%d" % i]
            A = relu(Z) * keep[:, [i]]
            cache["Z"].append(Z)
            cache["M"].append(keep[:, [i]])
            enc_out.append(A)
        H = np.concatenate(enc_out, axis=1)
        d1 = np.ones_like(H)
        if train and self.dropout > 0:
            d1 = (rng.rand(*H.shape) >= self.dropout) / (1 - self.dropout)
        Hd = H * d1
        Z2 = Hd @ P["V"] + P["v"]
        A2 = relu(Z2)
        d2 = np.ones_like(A2)
        if train and self.dropout > 0:
            d2 = (rng.rand(*A2.shape) >= self.dropout) / (1 - self.dropout)
        A2d = A2 * d2
        out = A2d @ P["U"] + P["u"]
        cache.update(H=H, d1=d1, Hd=Hd, Z2=Z2, d2=d2, A2d=A2d)
        return out, cache

    def backward(self, g, cache, wd=0.0):
        P = self.P
        G = {}
        G["U"] = cache["A2d"].T @ g + wd * P["U"]
        G["u"] = g.sum(0)
        dA2 = (g @ P["U"].T) * cache["d2"]
        dZ2 = dA2 * (cache["Z2"] > 0)
        G["V"] = cache["Hd"].T @ dZ2 + wd * P["V"]
        G["v"] = dZ2.sum(0)
        dH = (dZ2 @ P["V"].T) * cache["d1"]
        s = 0
        for i, b in enumerate(self.blocks):
            e = self.enc[i]
            dA = dH[:, s:s + e] * cache["M"][i]
            s += e
            dZ = dA * (cache["Z"][i] > 0)
            G["W%d" % i] = cache["X"][:, b].T @ dZ + wd * P["W%d" % i]
            G["c%d" % i] = dZ.sum(0)
        return G

    def predict(self, X, bs=4096):
        return np.vstack([self.forward(X[i:i + bs])[0]
                          for i in range(0, X.shape[0], bs)])


# ---------------- losses: (out, target) -> (loss, d loss / d out) -------
def loss_cox(out, tgt):
    t, e = tgt["t"], tgt["e"]
    lp = out[:, 0]
    ne = e.sum()
    g = np.zeros_like(out)
    if ne == 0:
        return 0.0, g
    m = lp.max()
    w = np.exp(lp - m)
    ut, inv = np.unique(t, return_inverse=True)
    wg = np.bincount(inv, weights=w, minlength=len(ut))
    S = np.cumsum(wg[::-1])[::-1]
    Si = S[inv]
    loss = -np.sum(e * (lp - m - np.log(Si))) / ne
    a = np.bincount(inv, weights=e / Si, minlength=len(ut))
    ca = np.cumsum(a)[inv]
    g[:, 0] = -(e - w * ca) / ne
    return loss, g


def sigmoid(z):
    return 0.5 * (1.0 + np.tanh(0.5 * z))


def loss_bce(out, tgt):
    Y, M = tgt["Y"], tgt["M"]
    p = np.clip(sigmoid(out), 1e-12, 1 - 1e-12)
    den = max(M.sum(), 1.0)
    loss = -np.sum(M * (Y * np.log(p) + (1 - Y) * np.log(1 - p))) / den
    return loss, (p - Y) * M / den


def hazard_targets(t, e, edges):
    """Interval targets for the logistic hazard: interval j covers
    (edges[j-1], edges[j]]. Y = 1 in the event's interval; M = 1 for the
    intervals the patient was at risk in."""
    n, K = len(t), len(edges)
    j = np.searchsorted(edges, t, side="left")
    j = np.minimum(j, K - 1)
    M = (np.arange(K)[None, :] <= j[:, None]).astype(float)
    Y = np.zeros((n, K))
    ev = (e == 1) & (t <= edges[-1])
    Y[np.where(ev)[0], j[ev]] = 1.0
    cens = ~ev
    last_full = (t[cens][:, None] >= edges[None, :])
    M[cens] = np.where(np.arange(K)[None, :] < j[cens][:, None], 1.0,
                       last_full.astype(float))
    return {"Y": Y, "M": M}


def loss_hazard(out, tgt):
    Y, M = tgt["Y"], tgt["M"]
    p = np.clip(sigmoid(out), 1e-12, 1 - 1e-12)
    n = out.shape[0]
    loss = -np.sum(M * (Y * np.log(p) + (1 - Y) * np.log(1 - p))) / n
    return loss, (p - Y) * M / n


def hazard_risk(out):
    h = sigmoid(out)
    return 1.0 - np.cumprod(1.0 - h, axis=1)


def softmax(z):
    z = z - z.max(1, keepdims=True)
    ez = np.exp(z)
    return ez / ez.sum(1, keepdims=True)


def deephit_targets(t, cause, edges, n_causes):
    """cause 0 = no event by its time; edges as in hazard_targets."""
    K = len(edges)
    j = np.minimum(np.searchsorted(edges, t, side="left"), K - 1)
    ev = (cause > 0) & (t <= edges[-1])
    return {"j": j, "c": np.where(ev, cause, 0), "t": t, "K": K,
            "C": n_causes}


def loss_deephit(out, tgt, alpha=0.5, sigma=0.1):
    alpha = tgt.get("alpha", alpha)
    sigma = tgt.get("sigma", sigma)
    K, C = tgt["K"], tgt["C"]
    j, c, t = tgt["j"], tgt["c"], tgt["t"]
    n = out.shape[0]
    P = softmax(out)
    cells = P[:, :C * K].reshape(n, C, K)
    F = np.cumsum(cells, axis=2)
    GP = np.zeros_like(P)
    GF = np.zeros((n, C, K))
    # negative log-likelihood
    ev = c > 0
    lik = np.empty(n)
    idx = np.where(ev)[0]
    lik[idx] = cells[idx, c[idx] - 1, j[idx]]
    cz = np.where(~ev)[0]
    surv = 1.0 - F[cz, :, :].sum(1)[np.arange(len(cz)), j[cz]]
    lik[cz] = surv
    lik = np.clip(lik, 1e-12, None)
    nll = -np.mean(np.log(lik))
    gcell = np.zeros((n, C, K))
    gcell[idx, c[idx] - 1, j[idx]] = -1.0 / lik[idx] / n
    for q, i in enumerate(cz):
        GF[i, :, j[i]] += 1.0 / lik[i] / n
    # ranking loss: an event of cause k at t_i should have a higher
    # cumulative incidence of k at t_i than anyone still event-free then
    rk = 0.0
    npair = 0
    GR = np.zeros((n, C, K))
    for k in range(1, C + 1):
        ii = np.where(c == k)[0]
        if not len(ii):
            continue
        Fi = F[ii, k - 1, j[ii]]
        Fo = F[:, k - 1, :][:, j[ii]].T
        valid = (t[None, :] > t[ii][:, None]).astype(float)
        eta = np.exp(-(Fi[:, None] - Fo) / sigma) * valid
        rk += eta.sum()
        npair += valid.sum()
        dFi = -(eta.sum(1)) / sigma
        dFo = eta / sigma
        np.add.at(GR, (ii, k - 1, j[ii]), dFi)
        for a, i in enumerate(ii):
            GR[:, k - 1, j[i]] += dFo[a]
    if npair > 0:
        rk = rk / npair
        GR = GR / npair
    loss = alpha * nll + (1 - alpha) * rk
    GF = alpha * GF + (1 - alpha) * GR
    # d/d cells: through the cumulative sums F (reverse cumulative)
    gcell = alpha * gcell + np.cumsum(GF[:, :, ::-1], axis=2)[:, :, ::-1]
    GP[:, :C * K] = gcell.reshape(n, C * K)
    g = P * (GP - (GP * P).sum(1, keepdims=True))
    return loss, g


def deephit_cif(out, n_causes, K):
    P = softmax(out)
    cells = P[:, :n_causes * K].reshape(len(P), n_causes, K)
    return np.cumsum(cells, axis=2)


LOSSES = {"cox": loss_cox, "bce": loss_bce, "hazard": loss_hazard,
          "deephit": loss_deephit}


def subset_target(tgt, idx):
    out = {}
    for k, v in tgt.items():
        out[k] = v[idx] if isinstance(v, np.ndarray) else v
    return out


# ---------------- training ----------------
def train(net, X, tgt, loss, lr=1e-3, wd=1e-4, bs=256, epochs=300,
          Xv=None, tv=None, patience=15, seed=0):
    """Adam on mini-batches. With validation data: early stopping on the
    validation loss, best weights kept; returns (best epoch, best val
    loss). Without: exactly `epochs` epochs."""
    rng = np.random.RandomState(seed)
    lf = LOSSES[loss]
    m = {k: np.zeros_like(v) for k, v in net.P.items()}
    v2 = {k: np.zeros_like(v) for k, v in net.P.items()}
    b1, b2, eps = 0.9, 0.999, 1e-8
    step = 0
    best, best_ep, best_P, wait = np.inf, 0, None, 0
    n = X.shape[0]
    for ep in range(1, epochs + 1):
        order = rng.permutation(n)
        for s in range(0, n, bs):
            bi = order[s:s + bs]
            if len(bi) < 2:
                continue
            out, cache = net.forward(X[bi], train=True, rng=rng)
            _, g = lf(out, subset_target(tgt, bi))
            G = net.backward(g, cache, wd)
            step += 1
            for k in net.P:
                m[k] = b1 * m[k] + (1 - b1) * G[k]
                v2[k] = b2 * v2[k] + (1 - b2) * G[k] ** 2
                mh = m[k] / (1 - b1 ** step)
                vh = v2[k] / (1 - b2 ** step)
                net.P[k] -= lr * mh / (np.sqrt(vh) + eps)
        if Xv is not None:
            vl = lf(net.predict(Xv), tv)[0]
            if vl < best - 1e-7:
                best, best_ep, wait = vl, ep, 0
                best_P = copy.deepcopy(net.P)
            else:
                wait += 1
                if wait >= patience:
                    break
    if Xv is not None and best_P is not None:
        net.P = best_P
        return best_ep, best
    return epochs, np.nan


def breslow_risk(lp_tr, t_tr, e_tr, lp_te, cuts):
    """Absolute risk 1 - exp(-H0(cut) exp(lp)) with the Breslow baseline
    of the training data (centred at the training maximum)."""
    o = np.argsort(t_tr)
    tt, ee, ll = t_tr[o], e_tr[o], lp_tr[o]
    w = np.exp(ll - ll.max())
    rs = np.cumsum(w[::-1])[::-1]
    out = np.zeros((len(lp_te), len(cuts)))
    for k, cu in enumerate(cuts):
        sel = (ee == 1) & (tt <= cu)
        if not sel.any():
            continue
        ut = np.unique(tt[sel])
        first = np.searchsorted(tt, ut, side="left")
        d = np.array([np.sum(sel & (tt == u)) for u in ut])
        H = np.sum(d / rs[first])
        out[:, k] = -np.expm1(-H * np.exp(lp_te - ll.max()))
    return out


def grad_check(seed=1):
    """Largest difference between analytic and numerical gradients for
    every loss on a small random problem, relative to the largest
    gradient checked (so near-zero entries cannot inflate it)."""
    rng = np.random.RandomState(seed)
    n, p = 24, 7
    X = rng.randn(n, p)
    t = np.round(rng.exponential(10, n)) + 0.5
    t[:3] = t[3]
    e = (rng.rand(n) < 0.6).astype(float)
    edges = np.array([3.5, 8.5, 15.5, 40.5])
    res = {}
    cases = {
        "cox": (1, {"t": t, "e": e}),
        "bce": (3, {"Y": (rng.rand(n, 3) < 0.4).astype(float),
                    "M": np.ones((n, 3))}),
        "hazard": (4, hazard_targets(t, e, edges)),
        "deephit": (2 * 4 + 1, deephit_targets(
            t, rng.randint(0, 3, n), edges, 2)),
    }
    for nm, (k, tgt) in cases.items():
        net = Net([np.arange(4), np.arange(4, p)], [5, 3], 6, k,
                  dropout=0.0, md=0.0, rng=np.random.RandomState(seed))
        out, cache = net.forward(X)
        _, g = LOSSES[nm](out, tgt)
        G = net.backward(g, cache, 0.0)
        diffs, mags = [], []
        for key in net.P:
            Pk = net.P[key]
            flat = Pk.reshape(-1)
            for q in rng.choice(len(flat), min(6, len(flat)), False):
                old = flat[q]
                flat[q] = old + 1e-6
                lp_ = LOSSES[nm](net.forward(X)[0], tgt)[0]
                flat[q] = old - 1e-6
                lm_ = LOSSES[nm](net.forward(X)[0], tgt)[0]
                flat[q] = old
                num = (lp_ - lm_) / 2e-6
                ana = G[key].reshape(-1)[q]
                diffs.append(abs(num - ana))
                mags.append(abs(ana))
        res[nm] = max(diffs) / max(1e-8, max(mags))
    return res