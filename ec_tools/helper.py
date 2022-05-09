r"""
Loose collection of helper functions used in several places in the package.
"""

import numpy as np

def find_x0_values(x, y, mode="all"):
    """Takes x and y as numpy arrays and returns the zero points by assuming linear behavior between
    the points.
    The mode determines if all or either positive or negative zero points should be returned.

    TODO:
    - add non-linear interpolation
    
    >>> x = np.array([10, 10.5, 11, 11.5, 12])
    >>> y = np.array([1, 1, -1, -1, 1])
    >>> find_x0_values(x,y)
    array([10.75, 11.75])
    """
    if mode == "all":
    	crossings = np.where(np.diff(np.sign(y)))[0]
    elif mode == "pos":
    	crossings = np.where(np.diff(np.sign(y)) > 0)[0]
    elif mode == "neg":
    	crossings = np.where(np.diff(np.sign(y)) < 0)[0]

    m = (y[crossings] - y[crossings+1]) / (x[crossings] - x[crossings+1])
    Δx = -y[crossings] / m 
    return Δx+x[crossings]