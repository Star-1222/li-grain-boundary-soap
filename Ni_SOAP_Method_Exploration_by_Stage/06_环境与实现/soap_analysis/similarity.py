"""SOAP kernels and similarity measures.

[NEW 2026-08-21] All functions assume row-normalized descriptor matrices
(N, D) as produced by soap_backend.soap_descriptors(..., normalize=True).
"""

from __future__ import annotations

import numpy as np


def _as_normalized(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"expected 2-D descriptor matrix, got {X.shape}")
    norms = np.linalg.norm(X, axis=1)
    norms[norms == 0.0] = 1.0
    return X / norms[:, None]


def kernel_matrix(X: np.ndarray, Y: np.ndarray = None) -> np.ndarray:
    """SOAP kernel (cosine) matrix K_ij = <x_i, y_j>."""
    Xn = _as_normalized(X)
    if Y is None:
        return Xn @ Xn.T
    Yn = _as_normalized(Y)
    return Xn @ Yn.T


def cosine_similarity(x: np.ndarray, y: np.ndarray) -> float:
    """Cosine similarity between two single descriptor vectors."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    nx = np.linalg.norm(x)
    ny = np.linalg.norm(y)
    if nx == 0.0 or ny == 0.0:
        return 0.0
    return float(np.dot(x, y) / (nx * ny))


def mean_vector(X: np.ndarray, normalize: bool = True) -> np.ndarray:
    """Mean SOAP vector of a structure, optionally re-normalized."""
    Xn = _as_normalized(X)
    m = Xn.mean(axis=0)
    if normalize:
        n = np.linalg.norm(m)
        if n > 0.0:
            m = m / n
    return m


def structure_similarity(Xa: np.ndarray, Xb: np.ndarray) -> float:
    """Similarity of two structures with MATCHED atoms (same N, same order).

    Returns the mean over atoms of the per-atom SOAP kernel
    (1.0 = identical environments, random ~ 0 for large descriptors).
    Use for repeated relaxations of the same initial model.
    """
    Xa = _as_normalized(Xa)
    Xb = _as_normalized(Xb)
    if Xa.shape != Xb.shape:
        raise ValueError(f"shape mismatch: {Xa.shape} vs {Xb.shape}")
    return float(np.mean(np.einsum("ij,ij->i", Xa, Xb)))


def fingerprint(X: np.ndarray) -> np.ndarray:
    """Rotation/permutation-insensitive structure fingerprint (mean SOAP)."""
    return mean_vector(X)
