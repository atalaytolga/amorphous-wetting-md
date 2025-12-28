# src/analysis/pbc.py
import numpy as np

def wrap_pbc(positions, box_edges):
    """
    Wraps coordinates into the centered box [-L/2, L/2].
    
    Parameters
    ----------
    positions : np.ndarray
        The coordinates to wrap.
    box_edges : np.ndarray
        The simulation box matrix (3x3)
        
    Returns
    -------
    np.ndarray
        Wrapped positions.
    """

    L = np.diag(box_edges)

    return ((positions + L/2) % L) - L/2