"""Chip-aware placement for QAOA on coupler-connected modular quantum devices."""
from .core import (
    cut_edges,
    exact_maxcut,
    exact_partition,
    kl_partition,
    make_graph,
    maxcut_diag,
    noisy_expectation,
    optimise_angles,
    qaoa_state,
    random_partition,
)

__all__ = [
    "maxcut_diag", "exact_maxcut", "qaoa_state", "optimise_angles", "noisy_expectation",
    "cut_edges", "random_partition", "kl_partition", "exact_partition", "make_graph",
]
