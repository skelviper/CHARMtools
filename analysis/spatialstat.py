"""
spatialstat.py  — 3D spatial (cross)correlation utilities

Reimplementation and correction of spatial autocorrelation in https://github.com/ammondongp/3D_ATAC_PALM

In this reimplementation, we fixed several issues in the original MATLAB code:
https://github.com/ammondongp/3D_ATAC_PALM/blob/92f0e15c275aef483e0a35dad75e7b2201d55e4a/spatialxcorr_3D_without_edge.m#L58


"""

from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from typing import Tuple, Optional, Callable

def shell_volume(r: float, dr: float) -> float:
    return 4.0 * np.pi * (r**2) * dr + (np.pi / 3.0) * (dr**3)

def _hull_of(points: np.ndarray) -> ConvexHull:
    return ConvexHull(points)

def _rejection_sample_in_hull(n: int,
                              hull: ConvexHull,
                              rng: np.random.Generator) -> np.ndarray:
    
    A = hull.equations[:, :3]
    b = hull.equations[:, 3]

    mins = np.min(hull.points, axis=0)
    maxs = np.max(hull.points, axis=0)
    box = maxs - mins
    out = []

    while len(out) < n:
        need = n - len(out)
        cand = rng.random((max(need * 2, 256), 3)) * box + mins
        inside = np.all((A @ cand.T + b[:, None]) <= 0, axis=0)
        if np.any(inside):
            chosen = cand[inside]
            out.extend(chosen[:need])
    return np.asarray(out, dtype=float)


# ==========================
# Ripley’s K/H (3D)
# ==========================

def ripleys_K_3D(points: np.ndarray,
                 r: float,
                 vol: Optional[float] = None) -> float:
    """
    3D Ripley's K。无向对 + n*(n-1) 
    K_Poisson(r) = (4/3) π r^3  
    """
    pts = np.asarray(points, dtype=float)
    n = len(pts)
    if n < 2:
        return 0.0

    if vol is None:
        hull = _hull_of(pts)
        vol = float(hull.volume)

    tree = cKDTree(pts)
    pairs = tree.query_pairs(r) 
    K = vol * (2.0 * len(pairs)) / (n * (n - 1))
    return float(K)


def ripleys_LH_3D(points: np.ndarray,
                  r: float,
                  vol: Optional[float] = None) -> Tuple[float, float]:
    """
    L(r) and H(r)：
    K_Poisson(r) = (4/3)π r^3
    L(r) = ( K / ((4/3)π) )^(1/3)
    H(r) = L(r) - r
    """
    K = ripleys_K_3D(points, r, vol=vol)
    L = (K / ((4.0 / 3.0) * np.pi)) ** (1.0 / 3.0)
    H = L - r
    return float(L), float(H)


# ==========================
# 3D Cross-correlation g_XY(r) with MC shape/edge correction
# ==========================

def _corr_shell_counts(X: np.ndarray,
                       Y: np.ndarray,
                       dr: float,
                       max_dist: float) -> Tuple[np.ndarray, np.ndarray]:

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    m = len(X)
    steps = int(np.floor(max_dist / dr))
    r = dr * (np.arange(steps) + 0.5)

    treeY = cKDTree(Y)
    counts = np.zeros(steps, dtype=np.int64)

    for i in range(steps):
        r_out = float(r[i] + dr / 2.0)
        r_in = float(max(0.0, r[i] - dr / 2.0))

        cnt_out = np.fromiter((len(treeY.query_ball_point(x, r_out)) for x in X),
                              dtype=np.int64, count=m)
        if r_in > 0.0:
            cnt_in = np.fromiter((len(treeY.query_ball_point(x, r_in)) for x in X),
                                 dtype=np.int64, count=m)
            cnt = cnt_out - cnt_in
        else:
            cnt = cnt_out
        counts[i] = int(cnt.sum())
    return counts, r


def spatialxcorr_3D_without_edge(X: np.ndarray,
                                 Y: np.ndarray,
                                 dr: float,
                                 vol: float,
                                 max_dist: float) -> Tuple[np.ndarray, np.ndarray, float]:

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    m = len(X)
    if m == 0 or len(Y) == 0:
        return np.zeros(0, dtype=float), np.zeros(0, dtype=float), 0.0

    rho_y = float(len(Y)) / float(vol)
    counts, r = _corr_shell_counts(X, Y, dr=dr, max_dist=max_dist)

    Corr = np.zeros_like(r, dtype=float)
    for i in range(len(r)):
        Vsh = shell_volume(float(r[i]), float(dr))
        Corr[i] = counts[i] / (m * rho_y * Vsh)
    return Corr, r, rho_y


def _align_zero(X: np.ndarray, Y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    mins = np.minimum(X.min(axis=0), Y.min(axis=0))
    return X - mins, Y - mins


def spatial_crosscorr_3d(X: np.ndarray,
                         Y: np.ndarray,
                         Nsim: int = 20,
                         dr: float = 0.5,
                         max_dist: float = 10.0,
                         vol: Optional[float] = None,
                         rng: Optional[np.random.Generator] = None
                         ) -> Tuple[np.ndarray, np.ndarray]:

    if rng is None:
        rng = np.random.default_rng()

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if len(X) == 0 or len(Y) == 0:
        return np.zeros(0, dtype=float), np.zeros(0, dtype=float)

    X, Y = _align_zero(X, Y)

    P = np.vstack([X, Y])
    hull = _hull_of(P)
    if vol is None:
        vol = float(hull.volume)

    Corr, r, _ = spatialxcorr_3D_without_edge(X, Y, dr=dr, vol=vol, max_dist=max_dist)

    CorrR = []
    for _ in range(Nsim):
        R1 = _rejection_sample_in_hull(len(X), hull, rng)
        R2 = _rejection_sample_in_hull(len(Y), hull, rng)
        corr_r, _, _ = spatialxcorr_3D_without_edge(R1, R2, dr=dr, vol=vol, max_dist=max_dist)
        CorrR.append(corr_r)
    Corr_R = np.mean(CorrR, axis=0)

    FCorr = Corr / Corr_R
    return FCorr, r


# ==========================
# 3D Autocorrelation g_XX(r) with MC baseline
# ==========================

def spatial_autocorr_3d(X: np.ndarray,
                        Nsim: int = 20,
                        dr: float = 0.5,
                        max_dist: float = 10.0,
                        vol: Optional[float] = None,
                        rng: Optional[np.random.Generator] = None
                        ) -> Tuple[np.ndarray, np.ndarray]:

    if rng is None:
        rng = np.random.default_rng()

    X = np.asarray(X, dtype=float)
    if len(X) == 0:
        return np.zeros(0, dtype=float), np.zeros(0, dtype=float)

    X0 = X.copy()
    FCorr, r = spatial_crosscorr_3d(X0, X0, Nsim=Nsim, dr=dr, max_dist=max_dist, vol=vol, rng=rng)
    return FCorr, r
