r"""
Implementation of a algorithm for semi-integration.
Fast Riemann-Liouville transformation (differintergration) - FRLT
based on
Pajkossy, T., Nyikos, L., 1984. Fast algorithm for differintegration. Journal of Electroanalytical Chemistry and Interfacial Electrochemistry 179, 65-69. https://doi.org/10.1016/S0022-0728(84)80275-2

"""

import numpy as np


def gruenwald(I, delta_x, v=-0.5):
    """
    Implementation of the Gruenwald algorithm for
    semi-integration (``v`` =-0.5) and semi-differentiation (``v`` =0.5)
    based on Oldham: Electrochemical Science and Technology, 2012

    EXAMPLES:

    Simple examples to compare the alg by a "double" semi-integration (i.e. resulting in a normal integration)
    with a numerical full integration from scipy. First one with linear graph for y:

    >>> from scipy.integrate import cumulative_trapezoid
    >>> x = np.linspace(0,1000, 1001)
    >>> delta_x = x[1] - x[0]
    >>> y = np.array([1]*1001)
    >>> np.allclose(gruenwald(gruenwald(y,delta_x),delta_x)[:-1], cumulative_trapezoid(y,x), rtol=1e-15)
    True

    Second test with more application-related values from a gaussian distribution function (from scipy):
    
    >>> from scipy.stats import norm
    >>> from scipy.integrate import cumulative_trapezoid
    >>> x = np.linspace(0, 8, 1001)
    >>> delta_x = x[1] - x[0]
    >>> y = norm.pdf(x,4,1)
    >>> np.allclose(gruenwald(gruenwald(y,delta_x),delta_x)[:-1], cumulative_trapezoid(y,x), rtol=1e-1)
    True

    """

    # No. of steps
    n_max = I.size
    # initialize with zeros
    g_1 = np.zeros(n_max)
    for N in range(1, n_max + 1):
        # value for n = N with w0 = 1
        g_i = I[0]
        #      go from N to 0
        for i in range(N - 1, 0, -1):
            g_i = g_i * (1 - (v + 1) / i) + I[N - i]
        g_1[N - 1] = g_i * np.sqrt(delta_x)
    return g_1

def riemann(y, delta_x, q=-0.5):
    """
    Implementation of the Riemann algorithm for
    semi-integration (``v`` =-0.5) and semi-differentiation (``v`` =0.5)
    based on
    Oldham: Electrochemical Science and Technology, 2012

    EXAMPLES:

    Simple examples to compare the alg by a "double" semi-integration (i.e. resulting in a normal integration)
    with a numerical full integration from scipy. First one with linear graph for y:

    >>> from scipy.integrate import cumulative_trapezoid
    >>> x = np.linspace(0,1000, 1001)
    >>> delta_x = x[1] - x[0]
    >>> y = np.array([1]*1001)
    >>> np.allclose(riemann(riemann(y,delta_x),delta_x)[:-1], cumulative_trapezoid(y,x), rtol=1e-0)
    True

    Second test with more application-related values from a gaussian distribution function (from scipy):

    >>> from scipy.stats import norm
    >>> from scipy.integrate import cumulative_trapezoid
    >>> x = np.linspace(0, 8, 1001)
    >>> delta_x = x[1] - x[0]
    >>> y = norm.pdf(x,4,1)
    >>> np.allclose(riemann(riemann(y,delta_x),delta_x)[:-1], cumulative_trapezoid(y,x), rtol=1e-0)
    True

    """

    # No. of steps
    n_max = y.size
    # initialize with zeros
    r_1 = np.zeros(n_max)

    if q == -0.5:
        sqrt_d_pi = (4 / 3) * np.sqrt(delta_x / np.pi)
    elif q == 0.5:
        sqrt_d_pi = 2 / np.sqrt(delta_x * np.pi)

    for N in range(1, n_max + 1):
        r_i = 0
        for i in range(1, N):
            r_i += y[i - 1] * (
                (N - i + 1) ** (1 - q) - 2 * (N - i) ** (1 - q) + (N - i - 1) ** (1 - q)
            )

        r_1[N - 1] = sqrt_d_pi * (
            y[N - 1]
            + y[0] * ((1 - q) * N ** (-q) - N ** (1 - q) + (N - 1) ** (1 - q))
            + r_i
        )

    return r_1


def fast_riemann(y, delta_x=1, q=-0.5, c1=8, c2=2):
    """
    Implementation of the fast Riemann algorithm for semi-integration.
    based on
    Pajkossy, T., Nyikos, L., 1984. Fast algorithm for differintegration. 
    Journal of Electroanalytical Chemistry and Interfacial Electrochemistry 179, 
    65-69. https://doi.org/10.1016/S0022-0728(84)80275-2
    
    Return the semiintegral R of order q for y with the x interval delta_x and the filter constants
    c1 and c2.

    Semi-integrating two times with order q = -0.5 should give the same result as integrating once.
    The relative error should not exceed 0.25 percent for 1000 and 0.5 percent per 10000 integration steps.
    
    EXAMPLES:

    >>> from scipy.integrate import cumulative_trapezoid
    >>> x = np.linspace(0,1000, 1001)
    >>> delta_x = x[1] - x[0]
    >>> y = np.array([1]*1001)
    >>> np.allclose(fast_riemann(fast_riemann(y, delta_x=delta_x), delta_x=delta_x), cumulative_trapezoid(y,x,initial=0), rtol=2.5e-03)
    True

    >>> x = np.linspace(0,1000, 10001)
    >>> delta_x = x[1] - x[0]
    >>> y = np.array([1]*10001)
    >>> np.allclose(fast_riemann(fast_riemann(y, delta_x=delta_x), delta_x=delta_x), cumulative_trapezoid(y,x,initial=0), rtol=5e-03)
    True
    
    """

    def prepare_kernel(q, delta_x, N, c1, c2):
        r"""
        Setup the integration kernel with the order q, the x interval delat_x, the length of the array N,
        and the filter constants c1 and c2.
        """
        tau0 = delta_x * N**0.5
        a0 = np.sin(np.pi * q) / (np.pi * q * tau0**q)
        # total number of filters
        n_filters = 2 * c1 * c2 + 1
        # dimension according to the number of filters
        # filter weights
        w1 = np.zeros(n_filters)
        w2 = np.zeros(n_filters)
        # auxiliary array
        s = np.zeros(n_filters)

        for i in range(2 * c1 * c2):
            j = i - c1 * c2
            a_j = (a0 / c2) * np.exp(j / c2)
            t_j = tau0 * np.exp(-j / (q * c2))
            w1[i] = t_j / (delta_x + t_j)
            w2[i] = a_j * (1 - w1[i])
        return s, w1, w2

    N = y.size
    R = np.zeros(N)
    s, w1, w2 = prepare_kernel(q, delta_x, N, c1, c2)
    for k in range(1, N):
        for i in range(s.size):
            s[i] = s[i] * w1[i] + y[k] * w2[i]
            R[k] = R[k] + s[i]
    return R
