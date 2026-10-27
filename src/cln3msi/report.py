"""Short per-phase summaries (markdown) and helpers for figures.

Every phase ends with summary.write(), which puts a markdown file in
<output_root>/summaries/<phase>_<mode>.md. Those are the files to paste back to me.
They list key numbers, warnings, and the figures written.
"""
from __future__ import annotations

import time

import matplotlib

matplotlib.use("Agg")   # no windows, works on a cluster node too
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402


def decimate(x, y, n: int = 20000):
    """Max-pool a long spectrum to ~n points for plotting (peaks stay visible, plots stay fast)."""
    import numpy as np
    if len(x) <= n:
        return x, y
    k = int(np.ceil(len(x) / n))
    m = len(x) // k * k
    yy = y[:m].reshape(-1, k)
    j = yy.argmax(axis=1) + np.arange(yy.shape[0]) * k
    return x[j], y[j]


class Summary:
    def __init__(self, ctx, phase: str, mode: str | None, title: str):
        self.ctx, self.phase, self.mode = ctx, phase, mode
        self.lines = [f"# {title}", "", f"_phase {phase}"
                      f"{', mode ' + mode if mode else ''}, written {time.strftime('%Y-%m-%d %H:%M')}_", ""]
        self.warnings: list[str] = []
        self.figures: list[str] = []

    def h(self, text: str):
        self.lines += ["", f"## {text}", ""]

    def p(self, text: str = ""):
        self.lines.append(text)

    def kv(self, d: dict):
        for k, v in d.items():
            if isinstance(v, float):
                v = f"{v:.4g}"
            self.lines.append(f"- **{k}**: {v}")

    def table(self, df: pd.DataFrame, max_rows: int = 40, floatfmt: str = ".7g"):
        if df is None or len(df) == 0:
            self.lines.append("_(empty table)_")
            return
        shown = df.head(max_rows).copy()
        # m/z with 4 decimals, p/q values in scientific notation, everything else compact
        for c in shown.columns:
            if not pd.api.types.is_float_dtype(shown[c]):
                continue
            name = str(c).lower()
            if "mz" in name:
                shown[c] = shown[c].map(lambda v: "" if pd.isna(v) else f"{v:.4f}")
            elif name.startswith(("p_", "q_", "p ", "q ")) or name in ("p", "q", "pval", "qval"):
                shown[c] = shown[c].map(lambda v: "" if pd.isna(v) else f"{v:.2e}")
            else:
                shown[c] = shown[c].map(lambda v: "" if pd.isna(v) else f"{v:.4g}")
        self.lines.append(shown.to_markdown(index=False, disable_numparse=True))
        if len(df) > max_rows:
            self.lines.append(f"\n_... {len(df) - max_rows} more rows in the CSV_")

    def warn(self, text: str, log=None):
        self.warnings.append(text)
        if log is not None:
            log.warning(text)

    def fig(self, fig, *parts: str, dpi: int = 130) -> str:
        p = self.ctx.out("figures", self.phase, *([self.mode] if self.mode else []), *parts)
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        rel = str(p.relative_to(self.ctx.output_root))
        self.figures.append(rel)
        return rel

    def write(self) -> str:
        body = list(self.lines)
        if self.warnings:
            body += ["", "## Warnings", ""] + [f"- {w}" for w in self.warnings]
        if self.figures:
            body += ["", "## Figures", ""] + [f"- `{f}`" for f in self.figures]
        name = f"{self.phase}{'_' + self.mode if self.mode else ''}.md"
        p = self.ctx.out("summaries", name)
        p.write_text("\n".join(body) + "\n")
        return str(p)
