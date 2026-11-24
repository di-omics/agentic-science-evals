"""Diffusion maps and diffusion pseudotime, written to follow the R package
destiny that the authors used (DiffusionMap, DPT, tips), plus small helpers.

destiny 2.x defaults reproduced here for a few hundred cells: every other cell
is a neighbour (k = n - 1), local kernel widths sigma_i = median distance to
the 5th, 6th and 7th nearest neighbours / 2, the Haghverdi et al. 2016 kernel
sqrt(2 s_i s_j / (s_i^2 + s_j^2)) exp(-d_ij^2 / (s_i^2 + s_j^2)), no self
transitions, density normalisation, and eigenvectors of the symmetric
conjugate rescaled to right eigenvectors of the transition matrix.

DPT(i, j)^2 = sum_k (lambda_k / (1 - lambda_k))^2 (psi_k(i) - psi_k(j))^2.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.distance import cdist


@dataclass
class DiffusionMap:
    evals: np.ndarray  # n_eigs non-trivial eigenvalues, descending
    psi: np.ndarray  # n x n_eigs right eigenvectors (diffusion components)
    sigma: np.ndarray  # local kernel widths
    data: np.ndarray  # the observations the map was built on

    def dpt(self, root: int) -> np.ndarray:
        w = self.evals / (1.0 - self.evals)
        diff = (self.psi - self.psi[root]) * w
        return np.sqrt((diff ** 2).sum(axis=1))

    def tips(self, start: int = 0) -> list[int]:
        """destiny's tip search: farthest from a start cell, farthest from that,
        then the cell farthest from both."""
        t1 = int(np.argmax(self.dpt(start)))
        d1 = self.dpt(t1)
        t2 = int(np.argmax(d1))
        d2 = self.dpt(t2)
        t3 = int(np.argmax(d1 + d2))
        return [t1, t2, t3]


def _local_sigma(dist: np.ndarray, n_local=(5, 6, 7)) -> np.ndarray:
    srt = np.sort(dist, axis=1)[:, 1:]  # drop self
    cols = [k - 1 for k in n_local]
    return np.median(srt[:, cols], axis=1) / 2.0


def diffusion_map(data: np.ndarray, n_eigs: int = 20) -> DiffusionMap:
    """Dense destiny-style diffusion map on the rows of `data`."""
    x = np.asarray(data, dtype=float)
    dist = cdist(x, x)
    sigma = _local_sigma(dist)
    s2 = sigma[:, None] ** 2 + sigma[None, :] ** 2
    w = np.sqrt(2.0 * np.outer(sigma, sigma) / s2) * np.exp(-(dist ** 2) / s2)
    np.fill_diagonal(w, 0.0)
    z = w.sum(axis=1)
    w = w / np.outer(z, z)  # density normalisation
    d = w.sum(axis=1)
    m = w / np.sqrt(np.outer(d, d))
    evals, evecs = np.linalg.eigh(m)
    order = np.argsort(-evals)[: n_eigs + 1]
    evals, evecs = evals[order], evecs[:, order]
    # fix the arbitrary sign of each eigenvector (largest-magnitude entry positive)
    signs = np.sign(evecs[np.argmax(np.abs(evecs), axis=0), np.arange(evecs.shape[1])])
    evecs = evecs * signs
    psi = evecs / np.sqrt(d)[:, None]
    return DiffusionMap(evals=evals[1:], psi=psi[:, 1:], sigma=sigma, data=x)


def spearman_dissimilarity(log_expr: np.ndarray) -> np.ndarray:
    """sqrt((1 - rho) / 2) between columns (cells) of a genes x cells matrix,
    as in the authors' R code (Spearman correlation between cells)."""
    ranks = np.apply_along_axis(_rank, 0, np.asarray(log_expr, dtype=float))
    ranks = ranks - ranks.mean(axis=0)
    ranks = ranks / np.sqrt((ranks ** 2).sum(axis=0))
    rho = np.clip(ranks.T @ ranks, -1.0, 1.0)
    return np.sqrt((1.0 - rho) / 2.0)


def _rank(v: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    return rankdata(v, method="average")


def spearman(a, b) -> float:
    from scipy.stats import spearmanr

    return float(spearmanr(np.asarray(a, float), np.asarray(b, float)).statistic)
