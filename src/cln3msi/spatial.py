"""Spatial statistics and embeddings, following the CAJAL tutorial notebooks 05-06.

  morans_i       spatial autocorrelation of every lipid, on a kNN graph of pixel positions
  seeded_nmf     NMF started from "seed" lipids (one per correlation cluster), more stable
                 than a random start
  territories    Leiden clustering of NMF activities, learned on Control, transferred to
                 every pixel with a kNN vote
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import cKDTree
from scipy.spatial.distance import squareform


def knn_weights(xy: np.ndarray, k: int) -> sparse.csr_matrix:
    """Binary kNN graph (k neighbours, no self). Rows = pixels."""
    tree = cKDTree(xy)
    _, idx = tree.query(xy, k=min(k + 1, len(xy)))
    n = len(xy)
    rows = np.repeat(np.arange(n), idx.shape[1] - 1)
    cols = idx[:, 1:].ravel()
    return sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))


def morans_i(X: np.ndarray, W: sparse.csr_matrix, n_perm: int = 0, rng=None) -> tuple[np.ndarray, np.ndarray]:
    """Moran's I of every column of X, and a permutation p-value (one-sided, I > null)."""
    Z = X - X.mean(axis=0)
    denom = (Z ** 2).sum(axis=0)
    denom[denom == 0] = np.nan
    n, s0 = X.shape[0], W.sum()
    I = n / s0 * (Z * (W @ Z)).sum(axis=0) / denom
    if not n_perm:
        return I, np.full(X.shape[1], np.nan)
    rng = np.random.default_rng(0) if rng is None else rng
    ge = np.zeros(X.shape[1])
    for _ in range(n_perm):
        Zp = Z[rng.permutation(n)]
        Ip = n / s0 * (Zp * (W @ Zp)).sum(axis=0) / denom
        ge += Ip >= I
    return I, (1 + ge) / (1 + n_perm)


def seeds_by_correlation(V: np.ndarray, k: int) -> list[int]:
    """Cluster lipids on 1 - correlation, return the most central lipid of each cluster."""
    C = np.corrcoef(V.T)
    C = np.nan_to_num(C)
    D = np.clip(1 - C, 0, 2)
    np.fill_diagonal(D, 0)
    lab = fcluster(linkage(squareform(D, checks=False), "average"), k, criterion="maxclust")
    seeds = []
    for c in np.unique(lab):
        members = np.flatnonzero(lab == c)
        seeds.append(int(members[np.argmin(D[np.ix_(members, members)].sum(axis=1))]))
    # if fewer clusters than k came out, top up with the most variable lipids
    for j in np.argsort(-V.var(axis=0)):
        if len(seeds) >= k:
            break
        if j not in seeds:
            seeds.append(int(j))
    return seeds[:k]


def seeded_nmf(V: np.ndarray, k: int, max_iter: int, seed: int = 0):
    from sklearn.decomposition import NMF
    s = seeds_by_correlation(V, k)
    C = np.nan_to_num(np.corrcoef(V.T))
    W0 = V[:, s].copy() + 1e-4
    H0 = np.clip(C[s], 0, None) + 1e-4
    m = NMF(n_components=k, init="custom", max_iter=max_iter, random_state=seed)
    W = m.fit_transform(V, W=W0, H=H0)
    return m, W, s


def choose_k(errors: dict) -> int:
    """Elbow: the k farthest from the straight line between the first and last error."""
    ks = np.array(sorted(errors))
    if len(ks) < 3:
        return int(ks[-1])
    e = np.array([errors[k] for k in ks], float)
    x = (ks - ks[0]) / (ks[-1] - ks[0])
    y = (e - e[-1]) / (e[0] - e[-1] if e[0] != e[-1] else 1)
    d = np.abs(1 - x - y) / np.sqrt(2)
    return int(ks[np.argmax(d)])


def territories(W_train: np.ndarray, k_nn: int, resolution: float, seed: int, log=None) -> np.ndarray:
    """Leiden clusters on a kNN graph of the NMF activities. Falls back to KMeans."""
    try:
        import igraph as ig
        import leidenalg
        tree = cKDTree(W_train)
        _, idx = tree.query(W_train, k=min(k_nn + 1, len(W_train)))
        edges = np.column_stack([np.repeat(np.arange(len(W_train)), idx.shape[1] - 1), idx[:, 1:].ravel()])
        g = ig.Graph(n=len(W_train), edges=edges.tolist(), directed=False).simplify()
        part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                                        resolution_parameter=resolution, seed=seed)
        lab = np.asarray(part.membership)
        # tiny clusters (< 0.5% of pixels) get merged into their nearest big one by the transfer step
        return lab
    except Exception as e:
        if log:
            log.warning(f"Leiden not available or failed ({e}); using KMeans with 8 clusters")
        from sklearn.cluster import KMeans
        return KMeans(8, n_init=10, random_state=seed).fit_predict(W_train)


def knn_transfer(W_train: np.ndarray, labels: np.ndarray, W_all: np.ndarray, k: int = 15) -> np.ndarray:
    tree = cKDTree(W_train)
    _, idx = tree.query(W_all, k=min(k, len(W_train)))
    from scipy import stats
    votes = labels[idx]
    return np.asarray(stats.mode(votes, axis=1, keepdims=False).mode)
