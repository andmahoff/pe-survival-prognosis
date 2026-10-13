"""suppress_small_cells.py
Prepares the results folder for publication.

Counts of 1 to 9 patients or events are replaced by "<10". Where the
rows of a table add up to a known total (deciles, sPESI classes, calendar
periods, subgroups) and only one cell in that set is small, the next
smallest count is replaced by "suppressed", so the small one cannot be
recovered by subtraction. Rates computed from a replaced count are
blanked, and shares that imply 1 to 9 patients are replaced by
"suppressed". A few text cells that quote small counts are reworded, and
Table 1 gets two more hidden cells (marked with a dagger) so its small
cell cannot be recovered from the row and column totals. Three tables
are left out because small counts could be recovered from them in many
rows. Absolute file paths are cut back to paths inside the project.
Every other cell is copied unchanged, and files without a rule are copied
byte for byte.

Usage:   python tools/suppress_small_cells.py <results in> <results out>
Outputs: the cleaned files in <results out>, and a list of the cells
         changed, with the reason but not the old value, in
         <results out>/../docs/suppressed_cells.csv
"""
import csv
import os
import re
import shutil
import sys

import openpyxl

SMALL = re.compile(r"^[1-9](\.0+)?$")

# left out of the public copy
WITHHELD = {
    "cal_curves.csv": "observed risk x n gives counts of 1 to 9 in many "
                      "bins; fin_cal_bins.csv merges bins to 10 or more "
                      "events",
    "s6_deciles.csv": "Kaplan-Meier risk x n gives counts of 1 to 9 in "
                      "several MIMIC-IV deciles",
    "lm_points.csv": "differences in the number at risk give the hidden "
                     "death counts",
}

# file: (count columns, columns that define a set of rows adding up to a
# known total, or None, rate columns blanked when a count is replaced)
COLUMN_RULES = {
    "conformal_groups.csv": (
        ["deaths_1y"],
        ["cohort", "alpha", "model", "method", "grouping"], []),
    "spesi_distribution.csv": (
        ["n", "deaths"], ["frame", "cohort", "horizon", "version"],
        ["death_rate"]),
    "subgroup_metrics.csv": (
        ["deaths"], ["frame", "model", "variable", "horizon", "measure"],
        []),
    "temporal_drift.csv": (["deaths_10"], ["cohort"], []),
    "lp_diag.csv": (["n_abs_z_gt5", "n_abs_z_gt10"], None, []),
    "lp_extreme_features.csv": (["n_extreme_patients", "n_patients"],
                                None, []),
}

# file: (n column, share columns); a share x n of 1 to 9 is replaced
SHARE_RULES = {
    "spesi_items.csv": ("n", ["age_over_80", "cancer",
                              "heart_failure_or_copd", "hr_110_or_more",
                              "sbp_under_100", "spo2_under_90"]),
    "stab_seed.csv": ("n", ["flip_at_10pct", "flip_at_20pct"]),
}

# long check files: (item column, value column, {item: new value}); a new
# value of None means "replace each count of 1 to 9 inside the text"
TEXT_RULES = {
    "audit_checks.csv": ("item", "value", {
        "dod ON discharge date, not in-hosp": "<10",
        "most index admissions, one patient": "<10",
        "deaths <=30 d: both/DISS only/native only": "suppressed",
        "deaths <=365 d: both/DISS only/native only": "suppressed"}),
    "bq_followup_checks.csv": ("item", "value", {
        "'DIED' location vs in-hospital death mismatch": "<10",
        "  location PSYCH FACILITY": "<10",
        "  location OTHER FACILITY": "<10",
        "  location ASSISTED LIVING": "<10",
        "Frame B, certainly in final year": "<10"}),
    "frames_summary.csv": ("item", "value", {
        "same-day CV and death within a year (CV)": "<10",
        "seed 42: one fold/patient; CV, death per fold": None}),
    "sens_checks.csv": ("item", "value", {
        "EHR+ECG: admissions kept / removed by the landmark":
            "3461 / 40 (deaths and CV readmissions; CV <10)",
        "3-mod: admissions kept / removed by the landmark":
            "1624 / 12 (deaths and CV readmissions; CV <10)"}),
}

# Table 1: the Home row is hidden in both 30-day columns, so the
# Rehabilitation count under 10 cannot be recovered from the totals
TABLE1_ROW = "Discharged to: Home, n (%)"
TABLE1_COLS = [2, 3]
TABLE1_NOTE_OLD = ("† marks a count not shown because it would reveal a "
                   "count under 10 in the other 30-day column.")
TABLE1_NOTE_NEW = ("† marks a count not shown because, with the row and "
                   "column totals, it would reveal a count under 10.")

WIN_PATH = re.compile(r"[A-Za-z]:\\[^,\"]*?\\SurvivalProject\\([^,\"]*)")


def read(path):
    raw = open(path, "rb").read()
    fmt = {"bom": raw.startswith(b"\xef\xbb\xbf"),
           "eol": "\r\n" if b"\r\n" in raw[:4096] else "\n"}
    with open(path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:], fmt


def write(path, header, rows, fmt=None):
    fmt = fmt or {"bom": False, "eol": "\n"}
    enc = "utf-8-sig" if fmt["bom"] else "utf-8"
    with open(path, "w", encoding=enc, newline="") as f:
        w = csv.writer(f, lineterminator=fmt["eol"])
        w.writerow(header)
        w.writerows(rows)


def is_small(v):
    return bool(SMALL.match(v.strip()))


def as_number(v):
    try:
        return float(v)
    except ValueError:
        return None


def column_rule(name, header, rows, log):
    counts, group_by, rates = COLUMN_RULES[name]
    idx = {h: i for i, h in enumerate(header)}
    hit = set()
    for col in counts:
        j = idx[col]
        if group_by is None:
            groups = {None: list(range(len(rows)))}
        else:
            groups = {}
            for r, row in enumerate(rows):
                key = tuple(row[idx[g]] for g in group_by)
                groups.setdefault(key, []).append(r)
        for members in groups.values():
            small = [r for r in members if is_small(rows[r][j])]
            if group_by is not None and len(small) == 1:
                rest = [r for r in members if r not in small and
                        (as_number(rows[r][j]) or 0) >= 10]
                if rest:
                    small.append(min(rest,
                                     key=lambda r: as_number(rows[r][j])))
            for r in small:
                if is_small(rows[r][j]):
                    why, mark = "count of 1 to 9", "<10"
                else:
                    why, mark = "complementary suppression", "suppressed"
                log.append((name, r + 2, col, why, mark))
                rows[r][j] = mark
                hit.add(r)
    for col in rates:
        j = idx[col]
        for r in sorted(hit):
            if rows[r][j] != "":
                log.append((name, r + 2, col, "rate from a replaced count",
                            "(blank)"))
                rows[r][j] = ""
    return rows


def share_rule(name, header, rows, log):
    n_col, shares = SHARE_RULES[name]
    jn = header.index(n_col)
    for col in shares:
        j = header.index(col)
        for r, row in enumerate(rows):
            n, v = as_number(row[jn]), as_number(row[j])
            if n is None or v is None:
                continue
            k = n * v
            if abs(k - round(k)) < 1e-6 and 1 <= round(k) <= 9:
                log.append((name, r + 2, col, "share giving 1 to 9 patients",
                            "suppressed"))
                row[j] = "suppressed"
    return rows


def text_rule(name, header, rows, log):
    item_col, value_col, items = TEXT_RULES[name]
    i, j = header.index(item_col), header.index(value_col)
    for r, row in enumerate(rows):
        if row[i] not in items:
            continue
        new = items[row[i]]
        if new is None:
            new = re.sub(r"(?<![\d.])[1-9](?![\d.])", "<10", row[j])
        if new != row[j]:
            log.append((name, r + 2, value_col, "count of 1 to 9 in text",
                        new))
            row[j] = new
    return rows


def table1_csv(name, header, rows, log):
    for r, row in enumerate(rows):
        if row[0] == TABLE1_ROW:
            for j in TABLE1_COLS:
                log.append((name, r + 2, header[j],
                            "complementary suppression", "†"))
                row[j] = "†"
    return rows


def table1_xlsx(src, dst, log):
    wb = openpyxl.load_workbook(src)
    ws = wb.active
    for row in ws.iter_rows():
        if row[0].value == TABLE1_ROW:
            for j in TABLE1_COLS:
                log.append(("table1.xlsx", row[j].row, row[j].coordinate,
                            "complementary suppression", "†"))
                row[j].value = "†"
        for c in row:
            if isinstance(c.value, str) and TABLE1_NOTE_OLD in c.value:
                c.value = c.value.replace(TABLE1_NOTE_OLD, TABLE1_NOTE_NEW)
                log.append(("table1.xlsx", c.row, c.coordinate,
                            "note on the dagger reworded", "(text)"))
    wb.save(dst)


def path_rule(name, header, rows, log):
    for r, row in enumerate(rows):
        for k, v in enumerate(row):
            if WIN_PATH.search(v):
                new = WIN_PATH.sub(
                    lambda x: x.group(1).replace("\\", "/"), v)
                log.append((name, r + 2, header[k], "absolute path", new))
                row[k] = new
    return rows


def main(src, dst):
    os.makedirs(dst, exist_ok=True)
    log = []
    for name in sorted(os.listdir(src)):
        p = os.path.join(src, name)
        if not os.path.isfile(p):
            continue
        if name in WITHHELD:
            log.append((name, "", "", "file left out: " + WITHHELD[name],
                        "(not published)"))
            continue
        if name == "table1.xlsx":
            table1_xlsx(p, os.path.join(dst, name), log)
            continue
        if not name.endswith(".csv"):
            shutil.copy2(p, os.path.join(dst, name))
            continue
        header, rows, fmt = read(p)
        n_before = len(log)
        if name in COLUMN_RULES:
            rows = column_rule(name, header, rows, log)
        if name in SHARE_RULES:
            rows = share_rule(name, header, rows, log)
        if name in TEXT_RULES:
            rows = text_rule(name, header, rows, log)
        if name == "table1.csv":
            rows = table1_csv(name, header, rows, log)
        rows = path_rule(name, header, rows, log)
        if len(log) > n_before:
            write(os.path.join(dst, name), header, rows, fmt)
        else:
            shutil.copy2(p, os.path.join(dst, name))
    out = os.path.join(os.path.dirname(os.path.abspath(dst)), "docs",
                       "suppressed_cells.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    write(out, ["file", "row", "column", "reason", "now"],
          [list(x) for x in log])
    print("changes:", len(log), "in", len({x[0] for x in log}), "files")
    print("list of changes:", out)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python tools/suppress_small_cells.py IN OUT")
    main(sys.argv[1], sys.argv[2])
