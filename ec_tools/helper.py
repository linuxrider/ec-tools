r"""
Loose collection of helper functions used in several places in the package.
"""

import numpy as np


def find_x0_values(x, y, mode="all"):
    """Takes x and y as numpy arrays and returns the ordered zero points by assuming linear behavior between
    the points.
    The mode determines if all or either positive or negative zero points should be returned.

    TODO:
    - add non-linear interpolation
    
    >>> x = np.array([10, 10.5, 11, 11.5, 12])
    >>> y = np.array([1, 1, -1, -1, 1])
    >>> find_x0_values(x, y)
    array([10.75, 11.75])

    >>> x = np.array([10, 10.5, 11, 11.5, 12])
    >>> y = np.array([1, 1, -1, -1, 1])
    >>> find_x0_values(x, y, mode='pos')
    array([11.75])

    >>> x = np.array([10, 10.5, 11, 11.5, 12])
    >>> y = np.array([1, 1, -1, -1, 1])
    >>> find_x0_values(x, y, mode='neg')
    array([10.75])

    Special case where y is exactly zero:
    >>> x = np.array([10, 10.5, 11, 11.5, 12])
    >>> y = np.array([1, 1, 0, -1, 1])
    >>> find_x0_values(x, y)
    array([11.  , 11.75])

    """
    signs = np.diff(np.sign(y))

    if mode == "all":
        exact_crossings = np.where((signs == 1) | (signs == -1))[0]
        non_exact_crossings = np.where((signs > 1) | (signs < -1))[0]
    elif mode == "pos":
        exact_crossings = np.where(signs == 1)[0]
        non_exact_crossings = np.where(signs > 1)[0]
    elif mode == "neg":
        exact_crossings = np.where(signs == -1)[0]
        non_exact_crossings = np.where(signs < -1)[0]

    m = (y[non_exact_crossings] - y[non_exact_crossings+1]) / (x[non_exact_crossings] - x[non_exact_crossings+1])
    Δx = -y[non_exact_crossings] / m

    return np.sort(np.concatenate([x[exact_crossings[1::2]], Δx + x[non_exact_crossings]]))

def determine_scan_rate(t, x):
    """Return scan rate of given t and x arrays.

    TEST:
    >>> t = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32])
    >>> E = np.array([0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.14, 0.13, 0.12, 0.11, 0.10, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14])
    >>> determine_scan_rate(t, E)
    0.005
    """

    return np.abs(np.diff(x) / np.diff(t)).mean()
