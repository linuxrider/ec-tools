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

    TEST:
    >>> y = np.array([1, 0.707106781186547, 0.577350269189626, 0.5, 0.447213595499958,
    ... 0.408248290463863, 0.377964473009227, 0.353553390593274, 0.333333333333333,
    ... 0.316227766016838, 0.301511344577764, 0.288675134594813, 0.277350098112615,
    ... 0.267261241912424, 0.258198889747161, 0.25, 0.242535625036333, 0.235702260395516,
    ... 0.229415733870562, 0.223606797749979, 0.218217890235992, 0.21320071635561,
    ... 0.208514414057075, 0.204124145231932, 0.2, 0.196116135138184, 0.192450089729875,
    ... 0.188982236504614, 0.185695338177052, 0.182574185835055, 0.179605302026775,
    ... 0.176776695296637, 0.174077655955698, 0.171498585142509, 0.169030850945703,
    ... 0.166666666666667, 0.164398987305357, 0.162221421130763, 0.160128153805087,
    ... 0.158113883008419, 0.156173761888606, 0.154303349962092, 0.152498570332605,
    ... 0.150755672288882, 0.149071198499986, 0.147441956154897, 0.145864991497895,
    ... 0.144337567297406, 0.142857142857143, 0.14142135623731
    ... ])
    >>> semi_integration(y)

    """
    N = y.size
    R = np.zeros(N)
    s, w1, w2 = prepare_kernel(q, Δx, N, c1, c2)
    for k in range(1, N):
        for i in range(s.size):
            s[i] = s[i]*w1[i] + y[k]*w2[i]
            R[k] = R[k] + s[i]
    return R

