"""SOAP analysis package for BCC Li grain boundaries.

[NEW 2026-08-21] Created for research item (e): verify that SOAP descriptors
can characterize equilibrated grain boundaries.  The package is
self-contained (numpy + scipy only) with optional dscribe / QUIP backends.

Modules
-------
soap_backend      SOAP power-spectrum descriptors (numpy reference +
                  optional dscribe / QUIP backends, periodic).
io_lammps         readers for LAMMPS data files and per-atom custom dumps.
bulk_reference    SOAP vector of the perfect BCC bulk at lattice constant a.
similarity        SOAP kernels, cosine similarity, structure fingerprints.
disorder_profile  per-atom structural disorder, x-binned GB profile, GB
                  width, disorder/energy correlation.
dedup             duplicate-structure detection for repeated quenches
                  (metastable-basin sampling).
compare_structures cross-GB comparison helpers (disorder vs theta).
run_analysis      command-line driver.
"""

__version__ = "1.0.0"

from .soap_backend import (BACKENDS, SoapHyper, available_backends,
                           power_spectrum_dims, soap_descriptors)
from .similarity import (cosine_similarity, kernel_matrix, mean_vector,
                         structure_similarity)
from .bulk_reference import bulk_bcc_positions, bulk_soap_reference
from .disorder_profile import (disorder_profile, gb_width_regions,
                               pearson_correlation, per_atom_disorder)

__all__ = [
    "BACKENDS", "SoapHyper", "available_backends", "power_spectrum_dims",
    "soap_descriptors", "cosine_similarity", "kernel_matrix", "mean_vector",
    "structure_similarity", "bulk_bcc_positions", "bulk_soap_reference",
    "disorder_profile", "gb_width_regions", "pearson_correlation",
    "per_atom_disorder", "__version__",
]
