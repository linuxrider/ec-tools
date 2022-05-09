r"""
Implementation of a algorithm for semi-integration.
Fast Riemann-Liouville transformation (differintergration) - FRLT
based on
Pajkossy, T., Nyikos, L., 1984. Fast algorithm for differintegration. Journal of Electroanalytical Chemistry and Interfacial Electrochemistry 179, 65–69. https://doi.org/10.1016/S0022-0728(84)80275-2
TODO: 
- evaluate if it is possible have varying Δx i.e. passing x array
"""

import numpy as np

def prepare_kernel(q, Δx, N, c1, c2):
    r"""
    Setup the integration kernel with the order q, the x interval Δx, the length of the array N,
    and the filter constants c1 and c2.
    """
    τ0 = Δx *N**0.5
    a0 = np.sin(np.pi * q) / (np.pi *q * τ0**q)
    # total number of filters
    n_filters = 2 * c1 * c2 + 1
    # dimension according to the number of filters
    # filter weights
    w1 = np.zeros(n_filters)
    w2 = np.zeros(n_filters)
    # auxiliary array
    s = np.zeros(n_filters)

    for i in range(2*c1*c2):
        j = i - c1 * c2
        a_j = (a0 / c2) * np.exp(j / c2)
        t_j = τ0 * np.exp(-j / (q*c2))
        w1[i] = t_j / (Δx + t_j)
        w2[i] = a_j * (1 - w1[i])
    return s, w1, w2

def semi_integration(y, q=-0.5, Δx=1, c1=8, c2=2):
    """
    Return the semiintegral R of order q for y with the x interval Δx and the filter constants
    c1 and c2.
    """
    N = y.size
    R = np.zeros(N)
    s, w1, w2 = prepare_kernel(q, Δx, N, c1, c2)
    for k in range(1, N):
        for i in range(s.size):
            s[i] = s[i]*w1[i] + y[k]*w2[i]
            R[k] = R[k] + s[i]
    return R

