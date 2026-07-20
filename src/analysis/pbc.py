import numpy as np


def wrap_pbc(positions, box_edges):
    """Wrap coordinates into the centered periodic box ``[-L/2, L/2)``.

    Parameters
    ----------
    positions : np.ndarray
        The coordinates to wrap.
    box_edges : np.ndarray
        The orthorhombic simulation box matrix.

    Returns
    -------
    np.ndarray
        Wrapped positions.
    """
    lengths = np.diag(box_edges)
    return ((positions + lengths / 2) % lengths) - lengths / 2
