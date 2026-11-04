"""Per-animal summaries and paired Mutant vs Control tests.

The unit of replication is the animal (one section each). So for every region (or
territory) we first collapse pixels to ONE number per animal and lipid: the mean normalized
intensity over that animal's pixels in the region. Only then do we test, with the 6 slides
as pairs (each slide = one Control + one Mutant).

  paired t-test on log2 values   primary. log2FC = mean of the per-slide differences,
                                 95% CI from the t distribution with n-1 df
  Wilcoxon signed-rank (exact)   sensitivity check. Heads-up: with 6 pairs the smallest
                                 possible two-sided p is 2/64 = 0.031, so it can't survive
                                 FDR across many lipids. It's there to spot t-tests driven
                                 by one weird pair, not to call hits.
  sign test (binomial)           counts pairs going up vs down, ignores the size of the change.
                                 Same floor as the Wilcoxon (0.031 with 6 pairs). Reported
                                 because "how many of the 6 pairs agree" is a natural question.
  sex difference (exploratory)   Welch t-test of the female vs male pair differences, i.e. the
                                 sex x genotype interaction. Sex itself cancels in the paired
                                 difference (both animals on a slide are the same sex). Sex is
                                 confounded with run here, so this is labelled exploratory.
  mixed model                    log2 ~ genotype + (1 | slide). Off by default
                                 (differential.run_mixedlm): with 2 animals per slide its
                                 variance estimate collapses and the p-values are far too small.

The paired t-test is exactly the fixed-block linear model log2 ~ genotype + slide, so
"paired t-test" and "linear model with slide as a block" are the same analysis.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def per_animal_means(X: np.ndarray, obs: pd.DataFrame, unit_mask: np.ndarray, min_pixels: int) -> tuple[pd.DataFrame, pd.Series]:
    """(animals x peaks) mean intensity over the pixels in unit_mask; NaN if too few pixels."""
    sub = obs.loc[unit_mask]
    n = sub.groupby("animal_id", observed=True).size()
    df = pd.DataFrame(X[unit_mask], index=sub["animal_id"].to_numpy()).groupby(level=0).mean()
    df.loc[n[n < min_pixels].index] = np.nan
    return df, n


def pseudocount(means: pd.DataFrame) -> np.ndarray:
    """Only needed for lipids where some animal's mean is exactly 0: then half the smallest
    positive mean. Otherwise 0, so the log2FC isn't squashed towards 0 by a big offset."""
    v = means.to_numpy()
    pos = np.where(v > 0, v, np.nan)
    with np.errstate(all="ignore"):
        c = np.nanmin(pos, axis=0) / 2
    has_zero = np.nanmin(np.where(np.isfinite(v), v, np.inf), axis=0) <= 0
    return np.where(has_zero & np.isfinite(c), c, 0.0)


def paired_tests(log2: pd.DataFrame, animals: pd.DataFrame, control: str, mutant: str,
                 run_mixed: bool = False, min_pairs: int = 2) -> pd.DataFrame:
    """log2: animals x peaks. animals: index animal_id with columns slide, genotype, sex."""
    info = animals.loc[log2.index]
    slides = sorted(info["slide"].unique())
    ctl_rows, mut_rows, sl = [], [], []
    for s in slides:
        a = info[(info.slide == s) & (info.genotype == control)].index
        b = info[(info.slide == s) & (info.genotype == mutant)].index
        if len(a) == 1 and len(b) == 1:
            ctl_rows.append(a[0]); mut_rows.append(b[0]); sl.append(s)  # noqa: E702
    D = log2.loc[mut_rows].to_numpy() - log2.loc[ctl_rows].to_numpy()     # pairs x peaks
    sex = info.loc[ctl_rows, "sex"].to_numpy()
    out = []
    for j, peak in enumerate(log2.columns):
        d = D[:, j]
        ok = np.isfinite(d)
        d = d[ok]
        n = len(d)
        rec = {"peak": peak, "n_pairs": n}
        if n >= max(2, min_pairs):
            m, sd = d.mean(), d.std(ddof=1)
            se = sd / np.sqrt(n)
            tcrit = stats.t.ppf(0.975, n - 1)
            rec.update(log2fc=m, ci_low=m - tcrit * se, ci_high=m + tcrit * se,
                       p_ttest=float(stats.ttest_1samp(d, 0).pvalue) if sd > 0 else (0.0 if m != 0 else 1.0),
                       n_pairs_up=int((d > 0).sum()))
            try:
                rec["p_wilcoxon"] = float(stats.wilcoxon(d, method="exact").pvalue) if np.any(d != 0) else 1.0
            except Exception:
                rec["p_wilcoxon"] = np.nan
            # sign test (binomial): only counts how many pairs go up. Same floor as Wilcoxon.
            nz = d[d != 0]
            rec["p_sign"] = float(stats.binomtest(int((nz > 0).sum()), len(nz)).pvalue) if len(nz) else 1.0
            # exploratory: is the effect different in females vs males? (confounded with run!)
            s_ok = sex[ok]
            groups = [d[s_ok == g] for g in np.unique(s_ok)]
            if len(groups) == 2 and all(len(g) >= 2 for g in groups):
                rec["p_sex_difference_exploratory"] = float(stats.ttest_ind(*groups, equal_var=False).pvalue)
                for g in np.unique(s_ok):
                    rec[f"log2fc_{g}"] = float(d[s_ok == g].mean())
        else:
            # too few pairs to test: report the effect, no p-value
            rec.update(log2fc=float(d.mean()) if n else np.nan, ci_low=np.nan, ci_high=np.nan,
                       p_ttest=np.nan, p_wilcoxon=np.nan, p_sign=np.nan, n_pairs_up=int((d > 0).sum()))
        if run_mixed and n >= max(3, min_pairs):
            rec["p_mixedlm"] = mixedlm_p(log2[peak], info, mutant)
        out.append(rec)
    res = pd.DataFrame(out)
    res["pairs"] = ",".join(sl)
    return res


def mixedlm_p(y: pd.Series, info: pd.DataFrame, mutant: str) -> float:
    import statsmodels.formula.api as smf
    df = pd.DataFrame({"y": y.to_numpy(), "mut": (info.loc[y.index, "genotype"] == mutant).astype(float).to_numpy(),
                       "slide": info.loc[y.index, "slide"].to_numpy()}).dropna()
    if df["mut"].nunique() < 2 or len(df) < 4:
        return np.nan
    # When the slide variance is ~0 the default optimizer hits a singular matrix; Nelder-Mead
    # still gets there. Note statsmodels reports a z-test here, which is a bit optimistic
    # with 12 animals, one more reason the paired t-test is the primary test.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for method in ("lbfgs", "nm"):
            try:
                fit = smf.mixedlm("y ~ mut", df, groups=df["slide"]).fit(reml=True, method=method)
                return float(fit.pvalues.get("mut", np.nan))
            except Exception:
                continue
    return np.nan


def bh(p: pd.Series) -> np.ndarray:
    p = p.to_numpy(float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q
