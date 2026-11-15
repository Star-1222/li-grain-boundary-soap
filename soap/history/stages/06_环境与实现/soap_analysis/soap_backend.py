"""SOAP power-spectrum descriptors with three interchangeable backends.

[NEW 2026-08-21]

Theory (Bartok et al., PRB 87, 184115 (2013)):

Around each atom i the neighbour density is expanded as

    rho_i(r) = sum_j exp(-|r - r_ij|^2 / (2 sigma^2))

and projected onto a radial basis {g_n} times spherical harmonics Y_lm,

    c^n_lm(i) = sum_j  <g_n Y_lm, exp(-|r - r_ij|^2 / (2 sigma^2))>

The (rotation-invariant) power spectrum is

    p^n n' l (i) = pi * sqrt(8 / (2l+1)) * sum_m c^n_lm(i) c^n'*_lm(i)

which we store for all unique pairs n <= n' and all l = 0..l_max, then
L2-normalize per atom.  The SOAP kernel between atoms is the dot product
of the normalized vectors (cosine similarity).

Backends
--------
numpy   Reference implementation of this module.  Gaussian radial basis
        with centers at rcut*(n+0.5)/n_max.  The atomic smearing sigma is
        folded into the radial projection analytically: projecting a
        Gaussian of width sigma onto a Gaussian basis function of width
        w_r gives a Gaussian in (r_j - c_n) of effective width
        sqrt(w_r^2 + sigma^2) (exact convolution of the two Gaussians).
        Requires only numpy + scipy.
dscribe The DScribe library (pip install dscribe), periodic SOAP with a
        GTO basis.  Same hyperparameters, but the basis conventions
        differ from the numpy reference, so descriptor VALUES are not
        comparable between backends -- only within one backend.
quippy  QUIP/libAtoms (quippy), the original SOAP implementation.
        Same caveat as dscribe.

backend="auto" picks the first available of dscribe -> quippy -> numpy.
For reproducibility pin the backend explicitly (e.g. --backend numpy).
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
from scipy.spatial import cKDTree

# spherical harmonics Y_lm(polar, azimuth).
# scipy >= 1.15: sph_harm_y(n, m, theta_polar, phi_azimuth)
# (physics convention; sph_harm_y(m, l, phi, theta) would SWAP the angles).
# older scipy: sph_harm(m, n, phi_azimuth, theta_polar).
try:  # scipy >= 1.15
    from scipy.special import sph_harm_y as _sph_harm_y

    def _ylm(l: int, m: int, theta: np.ndarray, phi: np.ndarray):
        """Y_lm with theta = polar angle, phi = azimuthal angle."""
        return _sph_harm_y(l, m, theta, phi)
except ImportError:  # pragma: no cover - older scipy
    from scipy.special import sph_harm as _sph_harm

    def _ylm(l: int, m: int, theta: np.ndarray, phi: np.ndarray):
        """Y_lm with theta = polar angle, phi = azimuthal angle."""
        return _sph_harm(m, l, phi, theta)


@dataclass(frozen=True)
class SoapHyper:
    """SOAP hyperparameters (all backends use the same names)."""
    rcut: float = 4.0        # cutoff radius [A]
    n_max: int = 6           # radial basis size
    l_max: int = 6           # maximum angular momentum
    sigma: float = 0.5       # atomic-density Gaussian width [A]


BACKENDS = ("auto", "dscribe", "quippy", "numpy")


def power_spectrum_dims(n_max: int, l_max: int) -> int:
    """Length of the per-atom SOAP vector: pairs(n<=n') x (l_max+1)."""
    return (n_max * (n_max + 1) // 2) * (l_max + 1)


# ----------------------------------------------------------------------------
# backend availability
# ----------------------------------------------------------------------------

def _try_import(module: str) -> Optional[object]:
    try:
        return __import__(module)
    except Exception:
        return None


def available_backends() -> List[str]:
    """Backends usable in the current environment (always includes numpy)."""
    out = ["numpy"]
    if _try_import("dscribe") is not None:
        out.insert(0, "dscribe")
    if _try_import("quippy") is not None:
        pos = out.index("dscribe") if "dscribe" in out else 0
        out.insert(pos, "quippy")
    return out


def _resolve_backend(backend: str) -> str:
    if backend not in BACKENDS:
        raise ValueError(f"unknown backend {backend!r}; choose from {BACKENDS}")
    if backend != "auto":
        return backend
    avail = available_backends()
    return avail[0]


# ----------------------------------------------------------------------------
# cell / input normalization
# ----------------------------------------------------------------------------

def _cell_matrix(cell) -> np.ndarray:
    cell = np.asarray(cell, dtype=float)
    if cell.shape == (3,):
        cell = np.diag(cell)
    if cell.shape != (3, 3):
        raise ValueError(f"cell must be 3-vector or 3x3, got shape {cell.shape}")
    return cell


def _symbols_list(symbols: Optional[Sequence[str]], n: int) -> List[str]:
    if symbols is None:
        return ["Li"] * n
    symbols = list(symbols)
    if len(symbols) == 1:
        return symbols * n
    if len(symbols) != n:
        raise ValueError(f"symbols length {len(symbols)} != number of atoms {n}")
    return symbols


# ----------------------------------------------------------------------------
# numpy reference backend
# ----------------------------------------------------------------------------

def _reference_soap(positions: np.ndarray, cell: np.ndarray, hp: SoapHyper,
                    normalize: bool) -> np.ndarray:
    n_atoms = positions.shape[0]
    n_dim = power_spectrum_dims(hp.n_max, hp.l_max)

    # --- periodic replica cloud -------------------------------------------
    # rows of `cell` are the lattice vectors.  For each dimension we need
    # ceil(rcut / L_i) replicas on each side so that every neighbour within
    # rcut is present in the cloud.
    lengths = np.array([np.linalg.norm(cell[i]) for i in range(3)])
    nrep = np.maximum(1, np.ceil(hp.rcut / lengths).astype(int))
    offsets = np.array(list(itertools.product(
        *[range(-int(n), int(n) + 1) for n in nrep])))          # (R, 3)
    shifts = offsets @ cell                                     # (R, 3)
    cloud = (positions[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
    tree = cKDTree(cloud)

    # --- angular / radial tables -------------------------------------------
    lm = [(l, m) for l in range(hp.l_max + 1) for m in range(-l, l + 1)]
    ls = np.array([l for l, _ in lm])
    radial_centers = hp.rcut * (np.arange(hp.n_max) + 0.5) / hp.n_max
    w_r = hp.rcut / (hp.n_max + 1)
    # exact convolution of the smearing Gaussian (sigma) with the radial
    # basis Gaussian (w_r): effective width sqrt(w_r^2 + sigma^2)
    w_eff = math.sqrt(w_r * w_r + hp.sigma * hp.sigma)
    pairs = [(n1, n2) for n1 in range(hp.n_max) for n2 in range(n1, hp.n_max)]

    X = np.zeros((n_atoms, n_dim))
    for i in range(n_atoms):
        idx = np.asarray(tree.query_ball_point(positions[i], hp.rcut), dtype=int)
        if idx.size == 0:
            continue
        dv = cloud[idx] - positions[i]
        r = np.linalg.norm(dv, axis=1)
        keep = r > 1e-12                      # drop the atom itself
        dv, r, idx = dv[keep], r[keep], idx[keep]
        if r.size == 0:
            continue
        theta = np.arccos(np.clip(dv[:, 2] / r, -1.0, 1.0))
        phi = np.arctan2(dv[:, 1], dv[:, 0])
        # smooth cutoff f_c(r) = 0.5 (cos(pi r / rcut) + 1)
        w = 0.5 * (np.cos(np.pi * r / hp.rcut) + 1.0)
        # radial projection (M, n_max)
        Rn = np.exp(-((r[:, None] - radial_centers[None, :]) ** 2)
                    / (2.0 * w_eff * w_eff)) * w[:, None]
        # spherical harmonics (M, n_lm)
        # [FIXED 2026-08-21] _ylm expects (l, m, theta_polar, phi_azimuth).
        # The old call _ylm(l, m, phi, theta) passed the two angles in the
        # wrong order, which broke the rotational invariance of the power
        # spectrum (rotated bulk similarity ~0.987 instead of ~1-1e-10).
        Y = np.empty((r.size, len(lm)), dtype=complex)
        for t, (l, m) in enumerate(lm):
            Y[:, t] = _ylm(l, m, theta, phi)
        C = Rn.T @ Y                           # (n_max, n_lm), complex
        # power spectrum over unique radial pairs
        vec = np.empty(n_dim)
        k = 0
        for l in range(hp.l_max + 1):
            sel = ls == l
            Cl = C[:, sel]
            P = np.real(Cl @ Cl.conj().T) * (math.pi * math.sqrt(8.0 / (2 * l + 1)))
            for (n1, n2) in pairs:
                vec[k] = P[n1, n2]
                k += 1
        X[i] = vec

    if normalize:
        norms = np.linalg.norm(X, axis=1)
        norms[norms == 0.0] = 1.0
        X = X / norms[:, None]
    return X


# ----------------------------------------------------------------------------
# dscribe backend
# ----------------------------------------------------------------------------

def _dscribe_soap(positions: np.ndarray, cell: np.ndarray, hp: SoapHyper,
                  symbols: List[str], normalize: bool) -> np.ndarray:
    from ase import Atoms
    from dscribe.descriptors import SOAP

    species = sorted(set(symbols))
    # [FIXED 2026-08-21] dscribe 2.x renamed the SOAP constructor keywords
    # (rcut -> r_cut, nmax -> n_max, lmax -> l_max) and REMOVED the
    # `normalize` option (the output is no longer L2-normalized, so we
    # normalize manually below).  Build the kwargs from the actual
    # signature so both dscribe 1.x and 2.x work.
    import inspect
    accepted = inspect.signature(SOAP.__init__).parameters
    kwargs = {}
    for new, old in (("r_cut", "rcut"), ("n_max", "nmax"), ("l_max", "lmax")):
        key = new if new in accepted else old
        if key in accepted:
            kwargs[key] = {"r_cut": hp.rcut, "n_max": hp.n_max,
                           "l_max": hp.l_max}[new]
    if "sigma" in accepted:
        kwargs["sigma"] = hp.sigma
    soap = SOAP(species=species, periodic=True, sparse=False, **kwargs)
    at = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    X = np.asarray(soap.create(at), dtype=float)
    if normalize:  # dscribe 2.x returns unnormalized vectors
        norms = np.linalg.norm(X, axis=1)
        norms[norms == 0.0] = 1.0
        X = X / norms[:, None]
    return X


# ----------------------------------------------------------------------------
# quippy backend
# ----------------------------------------------------------------------------

def _quippy_soap(positions: np.ndarray, cell: np.ndarray, hp: SoapHyper,
                 symbols: List[str], normalize: bool) -> np.ndarray:
    from ase import Atoms
    from quippy.descriptors import Descriptor

    species = sorted(set(symbols))
    species_str = " ".join(species)
    cmd = (f"soap cutoff={hp.rcut} n_max={hp.n_max} l_max={hp.l_max} "
           f"atom_sigma={hp.sigma} n_species={len(species)} "
           f"species_str='{species_str}'")
    if normalize:
        cmd += " normalise=T"
    desc = Descriptor(cmd)
    at = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    return np.asarray(desc.calc(at), dtype=float)


# ----------------------------------------------------------------------------
# public API
# ----------------------------------------------------------------------------

def soap_descriptors(positions, cell, symbols=None, *,
                     hyper: Optional[SoapHyper] = None,
                     backend: str = "auto",
                     normalize: bool = True) -> np.ndarray:
    """Per-atom SOAP power-spectrum descriptors of a periodic structure.

    Parameters
    ----------
    positions : (N, 3) array-like, atomic Cartesian positions [A].
    cell      : (3,) diagonal or (3, 3) array-like; rows are lattice vectors.
    symbols   : element symbols, one per atom (default: all "Li").
    hyper     : SoapHyper (default rcut=4, n_max=6, l_max=6, sigma=0.5).
    backend   : "auto" | "dscribe" | "quippy" | "numpy".
    normalize : L2-normalize each atom's descriptor (default True).

    Returns
    -------
    (N, power_spectrum_dims(n_max, l_max)) float array.
    """
    hp = hyper or SoapHyper()
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError(f"positions must be (N, 3), got {positions.shape}")
    cell = _cell_matrix(cell)
    symbols = _symbols_list(symbols, positions.shape[0])

    chosen = _resolve_backend(backend)
    if chosen == "numpy":
        return _reference_soap(positions, cell, hp, normalize)
    if chosen == "dscribe":
        return _dscribe_soap(positions, cell, hp, symbols, normalize)
    if chosen == "quippy":
        return _quippy_soap(positions, cell, hp, symbols, normalize)
    raise ValueError(f"unknown backend {chosen!r}")  # unreachable
