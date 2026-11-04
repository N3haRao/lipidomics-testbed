"""Phase 10: one self-contained HTML report (figures embedded, no internet needed to open it).

    python scripts/10_report.py --config config.yaml

Sections: data and QC, peaks and annotation, normalization, registration overlays,
lipid programs/territories, differential results (hypothesis classes first), side-by-side
ion images of the top hits (Control vs CLN3 KO per slide, shared colour scale), volcano
plots, caveats, and every phase summary at the end.

Writes <output_root>/report/report.html (+ the new ion-image figures under figures/10_report/).
"""
from __future__ import annotations

import base64
import html
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cln3msi import cli, provenance  # noqa: E402
from cln3msi.images import section_image  # noqa: E402
from cln3msi.report import Summary, plt  # noqa: E402
from cln3msi.util import read_h5ad  # noqa: E402

CSS = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1150px;margin:24px auto;padding:0 16px;
color:#1d1d1f;background:#fff;line-height:1.45}
h1{font-size:26px}h2{font-size:21px;border-bottom:1px solid #ddd;padding-bottom:4px;margin-top:36px}
h3{font-size:17px;margin-top:24px}
table{border-collapse:collapse;font-size:12px;margin:8px 0;display:block;overflow-x:auto}
td,th{border:1px solid #ddd;padding:3px 6px;text-align:right}th{background:#f4f4f6}
td:first-child,th:first-child{text-align:left}
img{max-width:100%;border:1px solid #eee;margin:6px 0}
.caveat{background:#fff6e5;border-left:4px solid #f0a020;padding:8px 12px;margin:8px 0}
details{margin:6px 0}summary{cursor:pointer;font-weight:600}
code{background:#f4f4f6;padding:1px 4px;border-radius:3px}
.small{color:#666;font-size:12px}
"""


def img_tag(path: Path, width: str = "100%") -> str:
    if not path.exists():
        return f"<p class='small'>(missing figure {html.escape(str(path.name))})</p>"
    b = base64.b64encode(path.read_bytes()).decode()
    return f"<img src='data:image/png;base64,{b}' style='width:{width}'>"


def table_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df is None or len(df) == 0:
        return "<p class='small'>(nothing to show)</p>"
    d = df.head(max_rows).copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            if str(c).startswith(("p_", "q_")):
                d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.2e}")
            elif "mz" in str(c):
                d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.4f}")
            else:
                d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.3g}")
    return d.to_html(index=False, escape=True, border=0)


def md_to_html(text: str) -> str:
    try:
        import markdown
        return markdown.markdown(text, extensions=["tables"])
    except Exception:
        return f"<pre>{html.escape(text)}</pre>"


def top_hits_images(ctx, mode: str, R: pd.DataFrame, S: Summary, n: int, log) -> Path | None:
    """Rows = hits, columns = slides x (Control, KO). Shared colour scale per hit."""
    p = ctx.out("processed", mode, "matrix_norm.h5ad", mkdir=False)
    if not p.exists() or R is None or R.empty:
        return None
    a = read_h5ad(p)
    obs = a.obs
    # best unit per peak, targeted hits first, then the rest by p
    # anatomical regions only (territories can exist on few slides), tested rows only
    best = R[(R["unit_kind"] == "region") & R["p_ttest"].notna()].sort_values("p_ttest").drop_duplicates("peak")
    hits = pd.concat([best[best["targeted"] & best["significant"]], best[~best["targeted"] & best["significant"]],
                      best[~best["significant"]]]).drop_duplicates("peak").head(n)
    if hits.empty:
        return None
    ctl, mut = ctx.cfg["project"]["genotypes"]["control"], ctx.cfg["project"]["genotypes"]["mutant"]
    slides = sorted(obs.slide.unique())
    fig, axes = plt.subplots(len(hits), 2 * len(slides), figsize=(1.25 * 2 * len(slides), 1.35 * len(hits) + 0.6),
                             squeeze=False)
    for r, (_, h) in enumerate(hits.iterrows()):
        j = a.var_names.get_loc(h["peak"])
        v = np.asarray(a.X[:, j]).ravel()
        vmax = np.percentile(v, 99.5)
        for k, s in enumerate(slides):
            for c, g in enumerate((ctl, mut)):
                ax = axes[r, 2 * k + c]
                sub = obs[(obs.slide == s) & (obs.genotype == g)]
                if len(sub):
                    img, *_ = section_image(obs, v, s, sub["section"].iloc[0])
                    ax.imshow(img, vmin=0, vmax=vmax, cmap=ctx.get("report.cmap", "viridis"), interpolation="nearest")
                ax.set_xticks([]); ax.set_yticks([])  # noqa: E702
                if r == 0:
                    ax.set_title(f"{s} {'Ctl' if g == ctl else 'KO'}", fontsize=7)
        lab = str(h.get("label", h["peak"]))
        axes[r, 0].set_ylabel(f"{lab[:24]}\n{h['unit']} log2FC {h['log2fc']:+.2f}", fontsize=6, rotation=0,
                              ha="right", va="center")
    fig.suptitle(f"{mode}: top hits, normalized intensity, one colour scale per row (Control | CLN3 KO per slide)",
                 fontsize=9)
    S.fig(fig, "top_hits_images.png", dpi=110)
    return ctx.output_root / S.figures[-1]


CAVEATS = [
    "<b>Genotype is confounded with position.</b> Control is always the top section and the Mutant always "
    "the bottom one. A technical top-to-bottom trend (matrix coating, instrument drift during the run) looks "
    "exactly like a genotype effect. We measured the trend in the off-tissue DAN signal (phase 03) and flag "
    "hits not clearly larger than it (<code>could_be_position</code>), but we cannot correct for it.",
    "<b>Sex is confounded with run.</b> All females were run in runs 1-3 and all males in run 4. Any sex "
    "comparison is exploratory only.",
    "<b>n = 6 pairs.</b> The animal is the unit (one section per animal); tests use per-animal means, paired by "
    "slide. With 6 pairs, power is limited and the Wilcoxon test can never reach p below 0.031.",
    "<b>Annotations are MS1 only, at sum-composition level.</b> No LC-MS/MS reference. Isobaric and isomeric "
    "species (e.g. BMP vs PG, GalCer vs GlcCer = HexCer, odd-chain vs ether lipids) can't be told apart; "
    "ambiguous matches are marked with *. Isotope (M+1, M+2) and type II overlaps are flagged.",
    "<b>No LC-MS validation</b> of abundances yet; MALDI intensities are semi-quantitative and ion "
    "suppression can differ between regions.",
    "<b>Registration is approximate.</b> The CCF has no true bregma; AP levels and region borders come from "
    "an automatic fit (check the overlays), and the KO brain may be atrophied.",
    "<b>The KO line is still to be confirmed</b> (whole-body vs conditional; IDs start with cKO). Outputs say "
    "'CLN3 KO' for now.",
]


def provenance_section(ctx) -> str:
    """Which code, config and packages produced these outputs (from metadata/provenance/)."""
    recs = [r for r in provenance.read_all(ctx) if r.get("phase") != "10_report"]
    if not recs:
        return ""
    T = pd.DataFrame([{"phase": r["phase"], "mode": r.get("mode") or "", "started": r.get("started"), "code_fingerprint": r.get("code_fingerprint"),
                       "git_commit": (r.get("git_commit") or "")[:12], "config_sha256": r.get("config_sha256"),
                       "sample_sheet_sha256": r.get("sample_sheet_sha256")} for r in recs]).sort_values(["phase", "mode"])
    now = provenance.collect(ctx)
    notes = []
    if T["code_fingerprint"].nunique() > 1:
        notes.append("The phases ran with <b>different code versions</b> (code_fingerprint differs). That is "
                     "fine if only later phases changed; rerun from the earliest changed phase to be sure.")
    if T["config_sha256"].nunique() > 1:
        notes.append("The phases ran with <b>different configs</b> (config_sha256 differs). Check that no setting "
                     "an earlier phase depends on was changed without rerunning it.")
    if not notes:
        notes.append("All phases ran with the same code and the same config.")
    pk = now.get("packages", {})
    pkg = ", ".join(f"{k} {v}" for k, v in pk.items() if v)
    return ("<h2>Reproducibility</h2><p class='small'>Every phase records the code fingerprint (a hash of src/, "
            "scripts/, reference/ and bin/, so it works for ZIP downloads too), the git commit if available, a hash "
            "of the merged config and of the sample sheet. The merged config of each phase is saved next to it in "
            "<code>metadata/provenance/</code>. To repeat this analysis, use code with the same fingerprint and the "
            "saved config.</p>"
            + "".join(f"<p>{n}</p>" for n in notes) + table_html(T, 20)
            + f"<p class='small'>Report written with pipeline v{html.escape(now['pipeline_version'])}, code "
              f"{html.escape(now['code_fingerprint'])}, Python {html.escape(now['python'])} on "
              f"{html.escape(now['platform'])}. Packages: {html.escape(pkg)}.</p>")


def slide_check_box(ctx, modes) -> str:
    """Top-of-report warning if 09 flagged a slide as a possible swap in any mode."""
    rows = []
    for mode in modes:
        p = ctx.output_root / "results" / "09_differential" / mode / "pair_consistency.csv"
        if p.exists():
            t = pd.read_csv(p)
            t.insert(0, "mode", mode)
            rows.append(t)
    if not rows:
        return ""
    T = pd.concat(rows, ignore_index=True)
    flagged = sorted(set(T.loc[T["flagged"].astype(bool), "slide"]))
    head = "<h2>Slide consistency check</h2>"
    txt = ("Each slide's log2(Mutant/Control) over all tested peaks is its difference profile. With correct "
           "labels the slides share the genotype effect, so profiles should correlate positively. Values are each "
           "slide's median Spearman correlation with the other slides.")
    if flagged:
        n_modes = {s: int(T.loc[T["flagged"].astype(bool) & (T["slide"] == s), "mode"].nunique()) for s in flagged}
        warn = ("<div class='caveat'><b>Possible label swap: slide(s) " + ", ".join(
            f"{html.escape(str(s))} ({n} of {len(rows)} modes)" for s, n in n_modes.items()) +
            ".</b> Their KO/Control difference is the mirror image of the other slides, which is what a swap "
            "of the two sections (sample sheet, or slide loaded the other way round) would look like. Results "
            "below use the sample sheet as it is. Check the sectioning/mounting records; only if they confirm "
            "the swap, list the slide in <code>qc.sections.flip_slides</code> and rerun from 03. The class-test "
            "CSVs have exploratory columns with the slide left out and flipped.</div>")
    else:
        warn = "<p>No slide flagged.</p>"
    return head + f"<p class='small'>{txt}</p>" + warn + table_html(T, 40)


def class_section(ctx, mode, fig) -> str:
    p = ctx.output_root / "results" / "09_differential" / mode / "class_tests.csv"
    if not p.exists():
        return ""
    R = pd.read_csv(p)
    pair_cols = [c for c in R.columns if c.startswith("pair_")]
    show = [c for c in ["unit", "feature", "kind", "n_species", "log2fc", "pct_change", "ci_low", "ci_high",
                        "p_ttest", "p_wilcoxon", "q_bh", "n_pairs_up"] + pair_cols if c in R]
    return ("<h3>Class totals and ratios (pre-specified)</h3><p class='small'>Species of a class added up per "
            "animal and region; ratios within the same animal and region. One BH family per mode "
            "(<code>q_bh</code>). pair_&lt;slide&gt; = that slide's log2(Mutant/Control).</p>"
            + table_html(R.sort_values("p_ttest")[show], 25)
            + img_tag(fig("09_differential", "class_tests_per_slide.png"), "85%")
            + img_tag(fig("09_differential", "class_tests_forest.png"), "70%")
            + img_tag(fig("09_differential", "pair_consistency.png"), "85%"))


def main():
    t0 = time.time()
    ctx, log, modes, _ = cli.start("10_report", __doc__)
    out = []
    P = ctx.cfg["project"]
    out.append(f"<h1>{html.escape(P['name'])}: MALDI-MSI lipidomics, {html.escape(P['mutant_label'])} vs Control</h1>")
    out.append(f"<p class='small'>Generated {time.strftime('%Y-%m-%d %H:%M')} from <code>{html.escape(str(ctx.config_path))}</code>."
               " Unpublished collaborator data, do not share.</p>")

    sheet = pd.read_csv(ctx.sample_sheet)
    out.append("<h2>Design</h2>")
    out.append(table_html(sheet[["animal_id", "genotype", "sex", "slide", "position_on_slide", "run_order"]], 20))
    out.append("<h2>Caveats</h2>" + "".join(f"<div class='caveat'>{c}</div>" for c in CAVEATS))
    out.append(slide_check_box(ctx, modes))
    out.append(provenance_section(ctx))

    files = ctx.output_root / "metadata" / "files.csv"
    if files.exists():
        f = pd.read_csv(files)
        cols = [c for c in ("mode", "file", "n_spectra", "grid_rows", "grid_cols", "storage_mode", "spectrum_type_data",
                            "polarity_header", "mz_min", "mz_max", "flags") if c in f]
        out.append("<h2>Data inventory</h2>" + table_html(f[cols], 30))

    for mode in modes:
        S = Summary(ctx, "10_report", mode, f"Report figures {mode}")
        fig = lambda ph, name: ctx.output_root / "figures" / ph / mode / name  # noqa: E731
        out.append(f"<h2>Mode {html.escape(mode)}</h2>")
        out.append("<h3>QC across slides</h3>")
        out.append(img_tag(fig("03_qc_all", "sections_all.png")))
        out.append(img_tag(fig("03_qc_all", "mass_error_all.png"), "60%"))
        out.append(img_tag(fig("03_qc_all", "gradient_all.png")))
        g = ctx.output_root / "results" / "03_qc_all" / mode / "gradient.csv"
        if g.exists():
            gd = pd.read_csv(g)
            out.append(table_html(gd[[c for c in gd.columns if c == "slide" or c.endswith(("pred_log2_bottom_vs_top", "spearman"))]]))
        out.append(img_tag(fig("03_qc_all", "mean_spectra_all.png")))
        out.append("<h3>Peaks and normalization</h3>")
        out.append(img_tag(fig("04_peaks_matrix", "peaks_on_mean_spectrum.png")))
        out.append(img_tag(fig("05_annotate", "annotation_overview.png")))
        out.append(img_tag(fig("06_normalize", "slide_effect_pca.png")))
        mc = ctx.output_root / "results" / "06_normalize" / mode / "method_comparison.csv"
        if mc.exists():
            out.append(table_html(pd.read_csv(mc)))
        out.append("<h3>Registration</h3>")
        out.append(img_tag(fig("07_register", "overlays.png")))
        rc = ctx.output_root / "results" / "07_register" / mode / "registration.csv"
        if rc.exists():
            out.append(table_html(pd.read_csv(rc)))
        out.append("<h3>Lipid programs and territories</h3>")
        out.append(img_tag(fig("08_explore", "programs_territories.png")))

        dpath = ctx.output_root / "results" / "09_differential" / mode / "differential.csv"
        R = pd.read_csv(dpath) if dpath.exists() else None
        out.append("<h3>Differential: hypothesis classes</h3>")
        if R is not None:
            show = [c for c in ("unit", "label", "hypothesis", "log2fc", "ci_low", "ci_high", "p_ttest", "q_targeted",
                                "p_wilcoxon", "n_pairs", "matches_hypothesis", "could_be_position", "ambiguous") if c in R]
            T = R[R["targeted"] & R["unit"].isin(["whole_section", "white_matter", "corpus_callosum", "wm_lipid",
                                                   "gm_lipid", "cortex", "hippocampus", "thalamus"]) & R["p_ttest"].notna()]
            out.append(table_html(T.sort_values("p_ttest")[show], 40))
            out.append(img_tag(fig("09_differential", "forest_hypothesis.png"), "75%"))
            out.append(class_section(ctx, mode, fig))
            out.append("<h3>Differential: untargeted, top per region</h3>")
            U = R[R.unit_kind == "region"].sort_values("p_ttest").groupby("unit").head(5)
            out.append(table_html(U[[c if c != "q_targeted" else "q_all" for c in show]], 60))
            out.append(img_tag(fig("09_differential", "volcano_regions.png")))
            ti = top_hits_images(ctx, mode, R, S, int(ctx.get("report.top_n_hits", 12)), log)
            if ti is not None:
                out.append("<h3>Top hits: Control vs CLN3 KO per slide</h3>" + img_tag(ti))
            n_sig = int(R["significant"].sum())
            n_pos = int(R.loc[R["significant"], "could_be_position"].fillna(False).astype(bool).sum())
            out.append(f"<div class='caveat'>{n_sig} significant tests in this mode, {n_pos} of them flagged "
                       "<code>could_be_position</code>.</div>")
        else:
            out.append("<p class='small'>(09 not run yet)</p>")
        S.write()

    out.append("<h2>All phase summaries</h2>")
    for p in sorted((ctx.output_root / "summaries").glob("*.md")):
        if p.name.startswith("10_report"):
            continue
        out.append(f"<details><summary>{html.escape(p.stem)}</summary>{md_to_html(p.read_text())}</details>")

    page = (f"<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,"
            f"initial-scale=1'><title>CLN3 MSI report</title><style>{CSS}</style></head><body>{''.join(out)}</body></html>")
    dest = ctx.out("report", "report.html")
    dest.write_text(page)
    log.info(f"wrote {dest} ({dest.stat().st_size / 1e6:.1f} MB)")
    cli.finish(log, t0)


if __name__ == "__main__":
    main()
