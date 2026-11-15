"""Per-atom structural disorder, GB profile and width, energy correlation.

[NEW 2026-08-21]

The "disorder" of atom i is its SOAP distance to the perfect-bulk
reference,

    d_i = 1 - <x_i, x_bulk>          in [0, 2]

d ~ 0 for bulk-like environments and is large at the grain boundary.
Binning d_i along x (the GB normal of the (hk0)[001] bicrystals) yields
the disorder profile, from which the GB width is estimated as the region
where the binned mean exceeds a fraction of the profile maximum.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np


def per_atom_disorder(X: np.ndarray, bulk_vector: np.ndarray) -> np.ndarray:
    """d_i = 1 - <x_i, x_bulk>; X and bulk_vector are L2-normalized."""
    X = np.asarray(X, dtype=float)
    if X.ndim != 2:
        raise ValueError(f"expected (N, D) matrix, got {X.shape}")
    b = np.asarray(bulk_vector, dtype=float).ravel()
    if b.shape[0] != X.shape[1]:
        raise ValueError(f"dimension mismatch: X {X.shape} vs bulk {b.shape}")
    nb = np.linalg.norm(b)
    if nb > 0.0:
        b = b / nb
    norms = np.linalg.norm(X, axis=1)
    norms[norms == 0.0] = 1.0
    Xn = X / norms[:, None]
    return 1.0 - Xn @ b


def disorder_profile(x_coords: np.ndarray, disorder: np.ndarray,
                     bin_width: float = 1.0, box_x: Optional[float] = None
                     ) -> Dict:
    """Bin the per-atom disorder along x.

    Parameters
    ----------
    x_coords  : (N,) atomic x coordinates [A].
    disorder  : (N,) per-atom disorder values.
    bin_width : bin width [A] (default 1.0).
    box_x     : periodic length in x (default: max(x) - min(x) + bin_width).

    Returns
    -------
    {"bin_centers", "mean_disorder", "max_disorder", "std_disorder",
     "counts", "bin_width", "profile_max"}
    """
    x = np.asarray(x_coords, dtype=float)
    d = np.asarray(disorder, dtype=float)
    if x.shape != d.shape:
        raise ValueError("x_coords and disorder must have the same length")
    if bin_width <= 0:
        raise ValueError("bin_width must be positive")

    if box_x is None:
        box_x = float(x.max() - x.min() + bin_width)
    n_bins = max(1, int(np.ceil(box_x / bin_width)))
    # periodic bin index (x may be wrapped or slightly negative)
    idx = np.floor(np.mod(x, box_x) / bin_width).astype(int)
    idx = np.clip(idx, 0, n_bins - 1)

    sums = np.bincount(idx, weights=d, minlength=n_bins)
    sqsums = np.bincount(idx, weights=d * d, minlength=n_bins)
    counts = np.bincount(idx, minlength=n_bins)

    mean = np.zeros(n_bins)
    std = np.zeros(n_bins)
    occupied = counts > 0
    mean[occupied] = sums[occupied] / counts[occupied]
    var = np.zeros(n_bins)
    var[occupied] = np.maximum(sqsums[occupied] / counts[occupied]
                               - mean[occupied] ** 2, 0.0)
    std[occupied] = np.sqrt(var[occupied])

    centers = (np.arange(n_bins) + 0.5) * bin_width
    return {"bin_centers": centers, "mean_disorder": mean,
            "std_disorder": std, "counts": counts,
            "bin_width": bin_width, "profile_max": float(mean.max())}


def gb_width_regions(profile: Dict, fraction: float = 0.5) -> List[Dict]:
    """Contiguous x-regions where mean_disorder > fraction * profile_max.

    Returns a list of {"x_lo", "x_hi", "width_A", "peak"} dicts, one per
    region (typically two, one per periodic GB plane).
    """
    mean = np.asarray(profile["mean_disorder"])
    centers = np.asarray(profile["bin_centers"])
    threshold = fraction * float(mean.max())
    above = mean > threshold

    regions: List[Dict] = []
    n = len(mean)
    # treat the bin array as periodic so a GB sitting at x ~ 0 is not split
    if above[0] and above[-1]:
        # rotate so the region containing bin 0 starts at index 0
        first_below = 0
        while first_below < n and above[first_below]:
            first_below += 1
        rot = first_below % n
        above = np.roll(above, -rot)
        centers = np.roll(centers, -rot)

    start: Optional[int] = None
    for i in range(n):
        if above[i] and start is None:
            start = i
        elif not above[i] and start is not None:
            regions.append((start, i - 1))
            start = None
    if start is not None:
        regions.append((start, n - 1))

    out: List[Dict] = []
    for lo, hi in regions:
        seg = mean[lo:hi + 1]
        xc = centers[lo:hi + 1]
        peak_i = int(np.argmax(seg))
        out.append({"x_lo": float(xc[0]), "x_hi": float(xc[-1]),
                    "width_A": float(xc[-1] - xc[0] + profile["bin_width"]),
                    "peak": float(seg[peak_i]),
                    "x_peak": float(xc[peak_i])})
    return out


def pearson_correlation(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson correlation coefficient; 0.0 if either array is constant."""
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if a.shape != b.shape:
        raise ValueError("arrays must have the same length")
    sa = a.std()
    sb = b.std()
    if sa == 0.0 or sb == 0.0:
        return 0.0
    return float(np.mean((a - a.mean()) * (b - b.mean())) / (sa * sb))
