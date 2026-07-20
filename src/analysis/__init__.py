"""
Analysis utilities for HALMD simulations.

This package contains analysis code used by scripts and notebooks.
"""

from .paths import DATA_DIR, FIGURES_DIR, PROCESSED_DIR, PROJECT_ROOT, RAW_DIR
from .pbc import wrap_pbc
from .plot_style import set_plot_style
from .simulation import Simulation

__all__ = [
    "PROJECT_ROOT",
    "DATA_DIR",
    "RAW_DIR",
    "PROCESSED_DIR",
    "FIGURES_DIR",
    "Simulation",
    "set_plot_style",
    "wrap_pbc",
]
