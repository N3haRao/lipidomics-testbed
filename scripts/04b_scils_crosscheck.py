"""Phase 04b (optional): compare our peak list with a SCiLS Lab export.

    python scripts/04b_scils_crosscheck.py --config config.yaml

Runs only if scils_crosscheck.enabled is true and a file is set for the mode. We don't
know how the exports are laid out, so this script first LOOKS at the file and prints what
it finds (separator, comment lines, columns, first rows), then tries two layouts:
  a) features in COLUMNS: column headers that parse as numbers in the mode's m/z range
  b) features in ROWS: a column called something like "m/z" with numbers in range
If neither works, it says so and shows you the preview, so you can tell me the layout.

Only the header and the first scils_crosscheck.preview_rows rows are read (plus the m/z
column in layout b), so a multi-GB export is fine.

Writes results/04b_scils_crosscheck/<mode>/{scils_preview.txt, matches.csv}, summary.
"""
from __future__ import annotations

import csv
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402

NUM = re.compile(r"[-+]?\d+(?:[.,]\d+)?")


def sniff_text(path: Path, n_lines: int = 60):
    raw = []
    with open(path, "r", errors="replace") as fh:
        for _ in range(n_lines):
            line = fh.readline()
            if not line:
                break
            raw.append(line.rstrip("\n"))
    comments = [ln for ln in raw if ln.startswith("#")]
    body = [ln for ln in raw if not ln.startswith("#") and ln.strip()]
    try:
        sep = csv.Sniffer().sniff("\n".join(body[:20]), delimiters=",;\t|").delimiter
    except Exception:
        sep = ";" if body and body[0].count(";") > body[0].count(",") else ","
    return raw, comments, sep


def to_float(s) -> float | None:
    m = NUM.fullmatch(str(s).strip()) or NUM.search(str(s))
    if not m:
        return None
    try:
        return float(m.group(0).replace(",", "."))
    except ValueError:
        return None


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("04b_scils_crosscheck", __doc__)
    C = ctx.cfg.get("scils_crosscheck", {})
    if not C.get("enabled"):
        log.info("scils_crosscheck.enabled is false: nothing to do")
        cli.finish(log, t0)
        return
    for mode in modes:
        f = (C.get("files") or {}).get(mode)
        if not f:
            log.info(f"{mode}: no SCiLS file configured, skipped")
            continue
        S = Summary(ctx, "04b_scils_crosscheck", mode, f"SCiLS export cross-check, {mode}")
        path = Path(f) if Path(f).is_absolute() else ctx.data_root / f
        if not path.exists():
            S.warn(f"file not found: {path}", log)
            S.write()
            continue
        n_prev = int(C.get("preview_rows", 20))
        S.h("What the file looks like")
        S.kv({"file": path.name, "size (MB)": round(path.stat().st_size / 1e6, 2)})
        preview_txt = []
        if path.suffix.lower() in (".xlsx", ".xls"):
            import openpyxl
            wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
            S.p(f"sheets: {wb.sheetnames} (reading the first one)")
            df = pd.read_excel(path, sheet_name=0, nrows=n_prev, header=None)
            wb.close()
            hdr_row = 0
            for i in range(min(10, len(df))):            # first row with >= 2 text cells = header
                if df.iloc[i].map(lambda v: isinstance(v, str)).sum() >= 2:
                    hdr_row = i
                    break
            df.columns = [str(c) for c in df.iloc[hdr_row]]
            df = df.iloc[hdr_row + 1:].reset_index(drop=True)
            reader = ("excel", hdr_row)
            preview_txt.append(df.head(n_prev).to_string())
        else:
            raw, comments, sep = sniff_text(path)
            S.p(f"separator: `{sep}` | comment lines at the top: {len(comments)}")
            preview_txt += raw[:30]
            df = pd.read_csv(path, sep=sep, comment=None, skiprows=len(comments), nrows=n_prev,
                             engine="python", on_bad_lines="skip")
            reader = ("csv", sep, len(comments))
        cols = [str(c) for c in df.columns]
        S.kv({"columns": len(cols), "first columns": ", ".join(cols[:12]) + (" ..." if len(cols) > 12 else "")})
        prev = ctx.out("results", "04b_scils_crosscheck", mode, "scils_preview.txt")
        prev.write_text("\n".join(preview_txt) + "\n\ncolumns:\n" + "\n".join(cols))
        log.info(f"preview written to {prev}")

        lo, hi = ctx.mode_def(mode).get("expected_mz", [0, 1e9])
        col_mz = [(c, to_float(c)) for c in cols]
        col_mz = [(c, v) for c, v in col_mz if v is not None and lo * 0.9 <= v <= hi * 1.1]
        scils_mz, layout = None, None
        if len(col_mz) >= 3:
            scils_mz = np.array([v for _, v in col_mz])
            layout = f"features in columns ({len(scils_mz)} m/z headers)"
        else:
            cand = [c for c in cols if re.search(r"m\s*/\s*z|\bmz\b|mass", c, re.I)]
            for c in cand:
                if reader[0] == "excel":
                    full = pd.read_excel(path, sheet_name=0, header=reader[1], usecols=[c])[c]
                else:
                    full = pd.read_csv(path, sep=reader[1], skiprows=reader[2], usecols=[c], engine="python")[c]
                v = pd.to_numeric(full.astype(str).str.replace(",", "."), errors="coerce").dropna().to_numpy()
                v = v[(v >= lo * 0.9) & (v <= hi * 1.1)]
                if len(v) >= 3:
                    scils_mz, layout = v, f"features in rows, column '{c}' ({len(v)} values)"
                    break
        S.h("Layout guess")
        if scils_mz is None:
            S.warn("couldn't find m/z values (neither numeric column headers nor an m/z column). "
                   f"Look at {prev.relative_to(ctx.output_root)} and tell me the layout.", log)
            S.p("```\n" + "\n".join(preview_txt[:15]) + "\n```")
            S.write()
            continue
        S.p(layout)

        mp = ctx.out("processed", mode, "peaks.csv", mkdir=False)
        if not mp.exists():
            S.warn("no peaks.csv from 04 yet; only the SCiLS side was read", log)
            S.write()
            continue
        ours = pd.read_csv(mp)["mz"].to_numpy()
        tol = float(C.get("match_ppm", 10))
        rows = []
        for m in np.sort(scils_mz):
            j = int(np.argmin(np.abs(ours - m)))
            e = (m - ours[j]) / ours[j] * 1e6
            rows.append({"scils_mz": m, "our_mz": ours[j], "ppm_diff": e, "matched": abs(e) <= tol})
        M = pd.DataFrame(rows)
        M.to_csv(ctx.out("results", "04b_scils_crosscheck", mode, "matches.csv"), index=False)
        ours_matched = {round(v, 6) for v in M.loc[M.matched, "our_mz"]}
        ours_only = [v for v in ours if round(v, 6) not in ours_matched]

        fig, ax = plt.subplots(figsize=(7, 3))
        ax.scatter(M["scils_mz"], M["ppm_diff"], c=np.where(M["matched"], "tab:green", "tab:red"), s=12)
        ax.axhspan(-tol, tol, color="0.9")
        ax.set_ylim(-5 * tol, 5 * tol)
        ax.set_xlabel("SCiLS m/z"); ax.set_ylabel("ppm to our nearest peak")  # noqa: E702
        S.fig(fig, "scils_vs_ours.png")

        S.h("Comparison")
        S.kv({"SCiLS features": len(M), f"matched within {tol:g} ppm": int(M.matched.sum()),
              "median ppm offset of matches (SCiLS - ours)": float(M.loc[M.matched, "ppm_diff"].median())
              if M.matched.any() else np.nan,
              "our peaks with no SCiLS feature": len(ours_only)})
        S.p("Our m/z values are recalibrated per slide in 04, so a constant offset vs SCiLS is expected "
            "if SCiLS wasn't recalibrated the same way.")
        S.table(M[~M.matched].head(30))
        log.info(f"summary: {S.write()}")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
