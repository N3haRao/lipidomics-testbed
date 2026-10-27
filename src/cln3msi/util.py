"""Odds and ends."""
from __future__ import annotations

import anndata as ad
import pandas as pd


def decat(df: pd.DataFrame) -> pd.DataFrame:
    """AnnData stores text columns as pandas Categoricals; turn them back into plain objects."""
    df = df.copy()
    for c in df.columns:
        if isinstance(df[c].dtype, pd.CategoricalDtype):
            df[c] = df[c].astype(object)
    return df


def read_h5ad(path) -> ad.AnnData:
    a = ad.read_h5ad(path)
    a.obs = decat(a.obs)
    a.var = decat(a.var)
    return a
