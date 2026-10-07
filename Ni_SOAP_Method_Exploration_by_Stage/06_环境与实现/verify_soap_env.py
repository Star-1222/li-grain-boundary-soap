#!/usr/bin/env python
"""verify_soap_env.py -- health check for the SOAP analysis environment.

[NEW 2026-08-21] Part of the 09_SOAP_Expansion_20260821 deliverable.

Checks, in order:
  1. Python version >= 3.9;
  2. required packages import (numpy, scipy);
  3. optional backends report their availability (dscribe, quippy);
  4. the scipy spherical-harmonics convention is handled correctly
     (guards the sph_harm -> sph_harm_y rename of scipy 1.15);
  5. a tiny SOAP computation reproduces the known bulk-BCC invariance:
     rotation about a GENERIC axis leaves the power spectrum unchanged
     (this fails if any angular convention is wrong);
  6. the SOAP package imports and reports its backends.

Exit code 0 = environment is good.

Usage:
    python environment/verify_soap_env.py
"""

from __future__ import annotations

import math
import sys

FAILURES = []


def check(name: str, ok: bool, detail: str = "") -> None:
    status = "OK  " if ok else "FAIL"
    print(f"[{status}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def main() -> int:
    print("=" * 72)
    print("SOAP environment verification (09_SOAP_Expansion_20260821)")
    print("=" * 72)

    # --- 1. python ---------------------------------------------------------
    v = sys.version_info
    check("Python >= 3.9", v >= (3, 9), f"{v.major}.{v.minor}.{v.micro}")

    # --- 2. core packages ----------------------------------------------------
    try:
        import numpy
        check("numpy", True, numpy.__version__)
    except Exception as exc:
        check("numpy", False, str(exc))
        # nothing else can run without numpy
        return _finish()

    try:
        import scipy
        import scipy.spatial  # noqa: F401  (cKDTree is required by the backend)
        check("scipy (+scipy.spatial.cKDTree)", True, scipy.__version__)
    except Exception as exc:
        check("scipy (+scipy.spatial.cKDTree)", False, str(exc))
        return _finish()

    # --- 3. optional backends -------------------------------------------------
    try:
        sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1]))
        from soap_analysis.soap_backend import available_backends
        backends = available_backends()
        check("SOAP backends available", True, ", ".join(backends))
    except Exception as exc:
        check("SOAP backends available", False, str(exc))
        return _finish()

    # --- 4. spherical harmonics convention ------------------------------------
    try:
        from soap_analysis.soap_backend import _ylm
        theta = [0.0, math.pi / 2, math.pi]
        phi = [1.234] * 3
        y10 = _ylm(1, 0, theta, phi)
        expected = [math.sqrt(3.0 / (4.0 * math.pi)) * math.cos(t) for t in theta]
        ok = all(abs(float(a.real) - b) < 1e-12 for a, b in zip(y10, expected))
        check("spherical harmonics convention Y_lm(polar, azimuth)", ok)
    except Exception as exc:
        check("spherical harmonics convention Y_lm(polar, azimuth)", False, str(exc))

    # --- 5. rotational invariance smoke test -----------------------------------
    try:
        import numpy as np
        from soap_analysis.bulk_reference import bulk_bcc_positions
        from soap_analysis.similarity import kernel_matrix
        from soap_analysis.soap_backend import SoapHyper, soap_descriptors

        hp = SoapHyper(rcut=4.0, n_max=4, l_max=4, sigma=0.5)   # small = fast
        pos, cell_vec = bulk_bcc_positions(a=3.5096, ncell=3)
        cell = np.diag(cell_vec)
        sym = ["Li"] * len(pos)

        axis = np.array([1.0, 2.0, 3.0])
        axis = axis / np.linalg.norm(axis)
        ang = math.radians(23.0)
        K = np.array([[0.0, -axis[2], axis[1]],
                      [axis[2], 0.0, -axis[0]],
                      [-axis[1], axis[0], 0.0]])
        R = np.eye(3) + math.sin(ang) * K + (1.0 - math.cos(ang)) * (K @ K)

        X = soap_descriptors(pos, cell, sym, hyper=hp, backend="numpy",
                             normalize=True)
        Xr = soap_descriptors(pos @ R.T, cell @ R.T, sym, hyper=hp,
                              backend="numpy", normalize=True)
        sim = float(np.diag(kernel_matrix(X, Xr)).min())
        check("SOAP rotational invariance (23 deg about [1,2,3])",
              sim > 1.0 - 1e-9, f"min cosine similarity = {sim:.12f}")
    except Exception as exc:
        check("SOAP rotational invariance (23 deg about [1,2,3])", False, str(exc))

    # --- 6. optional dscribe cross-check ----------------------------------------
    if "dscribe" in backends:
        try:
            from soap_analysis.soap_backend import soap_descriptors as sd
            Xd = sd(pos, cell, sym, hyper=hp, backend="dscribe",
                    normalize=True)
            ok = Xd.shape == X.shape
            check("dscribe backend computes full descriptor matrix", ok,
                  f"shape {Xd.shape}")
        except Exception as exc:
            check("dscribe backend computes full descriptor matrix", False,
                  str(exc))

    return _finish()


def _finish() -> int:
    print("-" * 72)
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} check(s) FAILED: {', '.join(FAILURES)}")
        return 1
    print("RESULT: environment is ready for SOAP analysis.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
