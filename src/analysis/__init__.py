"""
Analysis utilities for HALMD simulations.

This package contains analysis code used by scripts and notebooks.
"""

from .paths import PROJECT_ROOT, DATA_DIR, RAW_DIR, PROCESSED_DIR, FIGURES_DIR
from .pbc import wrap_pbc
from .simulation import Simulation

__all__ = [
    "PROJECT_ROOT",
    "DATA_DIR",
    "RAW_DIR",
    "PROCESSED_DIR",
    "FIGURES_DIR",
]
