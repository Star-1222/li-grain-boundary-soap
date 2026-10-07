"""SOAP reference vector of the perfect BCC bulk.

[NEW 2026-08-21]

Because the SOAP power spectrum is invariant under rotations, the bulk
reference does not depend on the crystal orientation -- only on the
lattice constant.  The per-atom structural disorder of a GB atom is then

    d_i = 1 - <x_i, x_bulk>

which is ~0 for perfect-bulk environments and grows with the local
structural distortion at the grain boundary.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .soap_backend import SoapHyper, soap_descriptors
from .similarity import mean_vector


def bulk_bcc_positions(a: float, ncell: int = 3) -> tuple:
    """Perfect BCC positions in an ncell^3 conventional periodic box.

    Returns (positions (2*ncell^3, 3), cell (3,)).
    """
    if a <= 0 or ncell < 2:
        raise ValueError(f"need a > 0 and ncell >= 2, got a={a} ncell={ncell}")
    pos = []
    for i in range(ncell):
        for j in range(ncell):
            for k in range(ncell):
                origin = np.array([i, j, k], dtype=float) * a
                pos.append(origin)
                pos.append(origin + np.array([0.5, 0.5, 0.5]) * a)
    pos = np.array(pos)
    cell = np.array([ncell * a, ncell * a, ncell * a])
    return pos, cell


def bulk_soap_reference(a: float, *, hyper: Optional[SoapHyper] = None,
                        backend: str = "auto",
                        element: str = "Li",
                        ncell: int = 3) -> np.ndarray:
    """L2-normalized mean SOAP vector of the perfect BCC bulk at lattice a."""
    pos, cell = bulk_bcc_positions(a, ncell=ncell)
    X = soap_descriptors(pos, cell, [element] * len(pos),
                         hyper=hyper, backend=backend, normalize=True)
    return mean_vector(X, normalize=True)
