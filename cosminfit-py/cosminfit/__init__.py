"""cosminfit: rate CFA fit against the COSMIN criterion.

Rates confirmatory factor analysis (CFA) results against the COSMIN
structural validity criterion "CFI or TLI or comparable measure > 0.95 OR
RMSEA < 0.06 OR SRMR < 0.08" (Prinsen et al. 2018; Mokkink et al. 2024)
and against alternative rules, and gives the population SRMR-RMSEA mapping.

This is a pure-Python port of the R package of the same name; ratings and
numeric results agree exactly with the R version (see the cross-check in
the source repository).

Main functions: :func:`rate_cfa`, :func:`rerate`, :func:`cosmin_rules`,
:func:`srmr_to_rmsea`, :func:`rmsea_to_srmr`, :func:`rmsea_equivalent_table`,
:func:`cfa_df`, :func:`n_moments`, :func:`wilson_ci`,
:func:`cluster_boot_prop`.
"""

from ._ci import cluster_boot_prop, wilson_ci
from ._mapping import (cfa_df, n_moments, rmsea_equivalent_table,
                       rmsea_to_srmr, srmr_to_rmsea)
from ._rerate import rerate
from ._rules import cosmin_rules, rate_cfa

__version__ = "0.1.0"

__all__ = [
    "cfa_df",
    "cluster_boot_prop",
    "cosmin_rules",
    "n_moments",
    "rate_cfa",
    "rerate",
    "rmsea_equivalent_table",
    "rmsea_to_srmr",
    "srmr_to_rmsea",
    "wilson_ci",
    "__version__",
]
