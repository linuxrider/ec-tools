r"""
The differentiation and the integration are common mathematical operations. 
The differentiation of an arbitrary function is often expressed by :math:`\frac{d}{dt} f(x)`.
Interestingly, the integration can also be expressed in a similar way:

:math:`\frac{d^{-1}}{dt^{-1}} f(x)= \int_0^t f(\tau) d\tau`

Therefore, these operations can be displayed in a more general way by:

:math:`\frac{d^{v}}{dt^{v}} f(x) = \int_0^t f(\tau) d\tau`

With:

* :math:`v= 1`: Differentiation

* :math:`v=-1`: Integration

Let's now introduce the so-called semi-operators. For :math:`v=1/2` we have the seimi differentiation and
(more interestingly for us) with :math:`v=-1/2` the seimi integration. 
The following picture visualize the idea of the semi integration and the semi differentiation.

.. image:: ../doc/files/images/semidif.png
  :width: 600
  :alt: Image 

The figure shows, that a semi integration of a peak function (bottom left) results in a hybrid function (top) and 
by performing another semi integration it brings a wave function (bottom right), which is equal to 
perform a "full" integration of the peak function. The opposide direction is similar, 
expect that a semi differentiation, respectively a "full" differentiation is performed instead.
  
**Semi Integration Methods**

Now we introduce some methods to apply the semi integration, resp. differentiation. 
These computations needs generally discrete values, i.e. the function graph (like above) has to be seperated into discrete finite values:

:math:`f(0), f(\delta), ..., f((N-1)\delta), f(N\delta)`

Here we assume that the step size :math:`\delta` is equidistant, meaning for a fixed set of x-Values N:

:math:`\delta = \frac{x_N}{N}`

The following algorithms (Gruenwald and Riemann & Liouville) are taken from Oldham in [1] and the fast Riemann from Pajkossy et. al. in [2].


Gruenwald Algorithms
^^^^^^^^^^^^^^^^^^^^
One sort of semi integration was introduced by Gruenwald [3] and 
Oldham shows in his web ressource 1244 from [1] how this type of semi integration can be applied as an algorithm, called G1. 

It can be generally expressed by taking the sum of the discrete function values multiplied with weights :math:`w_i` 
and then divided by the stepsize:

:math:`\frac{d^{\pm 0.5}}{dt^{\pm 0.5}} f(t) =\frac{1}{\delta^{\pm 0.5}} \sum_{n=0}^{N-1} w_n f(n\delta)`

The G1 algorithms is ideal for voltammograms like linear-scan or cyclic versions, where the early signals are small. 
**Note**, that these algorithms are less suitable for step and pulse varieties, in which the initial currents are large [1].

The weights can be expressed on different ways. The single weights :math:`w_i` often depends on their predecessor :math:`w_{i-1}`. 
As the factorial expression could lead to an overflow, this algorithm needs to be simplified.

The **Gruenwald G1 semi integration algorithm** is defined as follows:

:math:`\frac{d^{- 0.5}}{dt^{- 0.5}} f(t) \approx \sqrt{\delta} \sum_{n=1}^{N} w_{N-n} f(n\delta)`

Which can be also displayed in reverse summation to allow a more efficient
 implementation:

:math:`\frac{d^{- 0.5}}{dt^{- 0.5}} f(t) \approx \sqrt{\delta} \sum_{n=N}^{1} w_{N-n} f(n\delta)`

With:

* :math:`w_0 = 1`

* :math:`w_n = \frac{(n-0.5)w_{n-1}}{n} = (1-\frac{0.5}{n})w_{n-1}`

The previous definition can also be applied as  **Gruenwald G1 semi differentiation algorithm** with some changes:

:math:`\frac{d^{0.5}}{dt^{0.5}} f(t) \approx \frac{1}{\sqrt{\delta}} \sum_{n=0}^{N-1} w_{N-n} f(n\delta)`

With:

* :math:`w_0 = 1`

* :math:`w_n = \frac{(n-1.5)w_{n-1}}{n}`


Riemann and Liouville Algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another sort to determine the semi integral was introduced by Riemann and Liouville [4] and 
described by Oldham in [1] as R1 algorithm in his web ressource 1244. 
These sort of algorithms are mainly straightforward general-purpose algorithms.

The algorithm (from Web1242) for the **R1 semi integration** is:

:math:`\frac{d^{-1/2}}{dt^{-1/2}}f(t)=`

:math:`\frac{4}{3} \sqrt{\frac{\delta}{\pi}} \left[ f(N\delta) + \left\{ \frac{3}{2}\sqrt{N} - N^{3/2} + (N-1)^{3/2} \right\}f(0) + \right.`
:math:`\left.\sum_{n=1}^{N-1} \left\{ (N-n+1)^{3/2} - 2 (N-n)^{3/2} + (N-n-1)^{3/2} \right\}f(n\delta) \right]`

The R1 algorithms are not usable for application to currents that arise from potential steps or leaps [1], as:

1. The large current at t=0, immediately following the step is impossible to measure accurately and even if it would be possible,\
it is likely to be largely composed of a chemically uninteresting nonfaradaic component. But the algorithm still require a value of f(0).

2. The algorithm is based on the assumption that f(t) can be treated as an assemblage of linear segments, 
whereas faradaic currents arising from a potential step and are decidedly nonlinear with time. 

The general definition for the **R1 semi differentiation** is defined by:

:math:`\frac{d^{1/2}}{dt^{1/2}}f(t)=`

:math:`\frac{2}{\sqrt{\pi\delta}} \left[ f(N\delta) + \left\{ \frac{1}{2\sqrt{N}} - \sqrt{N} + \sqrt{N-1}\right\}f(0) +\right.`
:math:`\left. \sum_{n=1}^{N-1} \left\{ \sqrt{N-n+1} - 2 \sqrt{N-n} + \sqrt{N-n-1} \right\}f(n\delta) \right]`


Fast Riemann
^^^^^^^^^^^^

The following algorithm was introduced by Pajkossy et al. in 1984 [2] and is based on the Riemann-Liouville transformation (RLT). 
Its big advantage is, that the computation time increases only linearly with the number of points (N). 
Here it is necessary to define some input variables (beside the t & I(t) data), where q is equal to v, describing a semi integration or semi differentiation. 
:math:`\Delta_t` defines the constant time intervall (i.e. :math:`t_2 - t_1`) and :math:`c_1, c_2` are constant values,
which we set by default to :math:`c_1=8, c_2=2`, as Pajkossy recommends (but are still changeable). With these variables, the procedure of the algorithm can be described as pseudo code:

**Input** 
:math:`q, N, \Delta_t, c_1, c_2, I`

:math:`t_0 = \Delta_t N^{1/2}` 

:math:`a_0 = \sin(\pi q)/(\pi qt_0^q)`

**For** :math:`i = 0,2c_1c_2`

    :math:`j = i-c_1c_2`
    
    :math:`a_j = (a_0/c_2)\exp(j/c_2)`

    :math:`t_j = t_0 \exp(-j/qc_2)`

    :math:`w_1(i) = t_j / (\Delta_t + t_j)`

    :math:`w_2(i) = a_j(1-w_1(i))`

    :math:`s(i) = 0`

**For** :math:`k=1,N`

    :math:`R(k)=0`
    
    **For** :math:`i=0,2c_1c_2`
    
        :math:`s(i)= s(i)w_1(i) + I(k)w_2(i)`

        :math:`R(k) = R(k) + s(i)`

Here, :math:`R` represents the calculated semi integral, i.e. :math:`R \approx \frac{d^{v}}{dt^{v}} I(t)`.

"""
import numpy as np
from transonic import jit


def semi_integration(y, x, v=-0.5, alg="frlt", transonic_backend="pythran", d_tol=1e-5):
    r"""
    To perform a semi integration (``v`` =-0.5) or semi differentiation (``v`` =0.5)
    on given data (``y`` and ``t`` ) with speed up by transonic (with numba or pythran backend) or without simply using python
    different methods are implemented.

    Available algorithms (``alg`` ):

    ``frlt``: Fast Riemann-Liouville transformation (default)

    ``g1``: Gruenwald

    ``r1``: Riemann and Liouville

    Available settings (``transonic_backend`` ):

    ``python``: Transonic package with python backend (default)

    ``numba``: Transonic package with numba backend

    ``pythran``: Transonic package with pythran backend

    If steps are not equally spaced, ``d_tol`` (by default: 1e-5) defines the maximum step size difference on average.
    """

    # Calc average step size
    deltas = np.diff(x)
    delta_x = deltas.mean()

    # warning if avg. delta_x differs to much from single ones
    if np.max(np.abs(deltas / delta_x - 1)) >= d_tol:
        print("Warning: step size tolerance reached!")

    if alg not in ["g1", "r1", "frlt"] or transonic_backend not in [
        "python",
        "numba",
        "pythran",
    ]:
        raise ValueError("\nNo matching setting for alg or transonic_backend")

    if alg == "frlt":
        if transonic_backend == "python":

            @jit(backend="python")
            def fast_riemann_python(y, delta_x, v):
                return fast_riemann(y, delta_x, v)

            return fast_riemann_python(y, delta_x, v)

        if transonic_backend == "numba":

            @jit(backend="numba")
            def fast_riemann_numba(y, delta_x, v):
                return fast_riemann(y, delta_x, v)

            return fast_riemann_numba(y, delta_x, v)

        if transonic_backend == "pythran":

            @jit(backend="pythran")
            def fast_riemann_pythran(y, delta_x, v):
                return fast_riemann(y, delta_x, v)

            return fast_riemann_pythran(y, delta_x, v)

    elif alg == "g1":
        if transonic_backend == "python":

            @jit(backend="python")
            def gruenwald_python(y, delta_x, v):
                return gruenwald(y, delta_x, v)

            return gruenwald_python(y, delta_x, v)

        if transonic_backend == "numba":

            @jit(backend="numba")
            def gruenwald_numba(y, delta_x, v):
                return gruenwald(y, delta_x, v)

            return gruenwald_numba(y, delta_x, v)

        if transonic_backend == "pythran":

            @jit(backend="pythran")
            def gruenwald_pythran(y, delta_x, v):
                return gruenwald(y, delta_x, v)

            return gruenwald_pythran(y, delta_x, v)

    elif alg == "r1":
        if transonic_backend == "python":

            @jit(backend="python")
            def riemann_python(y, delta_x, v):
                return riemann(y, delta_x, v)

            return riemann_python(y, delta_x, v)

        if transonic_backend == "numba":

            @jit(backend="numba")
            def riemann_numba(y, delta_x, v):
                return riemann(y, delta_x, v)

            return riemann_numba(y, delta_x, v)

        if transonic_backend == "pythran":

            @jit(backend="pythran")
            def riemann_pythran(y, delta_x, v):
                return riemann(y, delta_x, v)

            return riemann_pythran(y, delta_x, v)


def gruenwald(y, delta_x, v=-0.5):
    """
    Gruenwald Algorithm
    ^^^^^^^^^^^^^^^^^^^

    Implementation of the Gruenwald algorithm for
    semi-integration (``v`` =-0.5) and semi-differentiation (``v`` =0.5)
    based on Oldham: [1]

    Input:

    ``y``: y-values

    ``delta_x``: step size (i.e. x2-x1)

    ``v``: -0.5 (default) or in range -1 < v < 1

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
    if v != (-0.5 or 0.5):
        print("\nWarning: algorithm is only tested for v=0.5 and v=-0.5.")
        print("         Other values for v are not verified!\n")

    # No. of steps
    n_max = y.size
    # initialize with zeros
    g_1 = np.zeros(n_max)
    for N in range(1, n_max + 1):
        # value for n = N with w0 = 1
        g_i = y[0]
        #      go from N to 0
        for i in range(N - 1, 0, -1):
            g_i = g_i * (1 - (v + 1) / i) + y[N - i]
        g_1[N - 1] = g_i * np.sqrt(delta_x)
    return g_1


def riemann(y, delta_x, v=-0.5):
    """

    Riemann Algorithm
    ^^^^^^^^^^^^^^^^^

    Implementation of the Riemann algorithm for
    semi-integration (``v`` =-0.5) and semi-differentiation (``v`` =0.5)
    based on Oldham: [1]

    Input:

    ``y``: y-values

    ``delta_x``: step size (i.e. x2-x1)

    ``v``: -0.5 (default) or 0.5

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

    if v not in (0.5, -0.5):
        raise ValueError(
            "\nError: This algorithm accept right now only v=0.5 and v=-0.5."
        )

    # No. of steps
    n_max = y.size
    # initialize with zeros
    r_1 = np.zeros(n_max)

    if v == -0.5:
        sqrt_d_pi = (4 / 3) * np.sqrt(delta_x / np.pi)
    elif v == 0.5:
        sqrt_d_pi = 2 / np.sqrt(delta_x * np.pi)

    for N in range(1, n_max + 1):
        r_i = 0
        for i in range(1, N):
            r_i += y[i - 1] * (
                (N - i + 1) ** (1 - v) - 2 * (N - i) ** (1 - v) + (N - i - 1) ** (1 - v)
            )

        r_1[N - 1] = sqrt_d_pi * (
            y[N - 1]
            + y[0] * ((1 - v) * N ** (-v) - N ** (1 - v) + (N - 1) ** (1 - v))
            + r_i
        )

    return r_1


def fast_riemann(y, delta_x=1, q=-0.5, c1=8, c2=2):
    """

    Fast Riemann Algorithm
    ^^^^^^^^^^^^^^^^^^^^^^

    Implementation of the fast Riemann algorithm for semi-integration.
    based on Pajkossy et al [2]. Return the semiintegral R of order q for y
    with the x interval delta_x and the filter constants c1 and c2.

    Input:

    ``y``: y-values

    ``delta_x``: step size (i.e. x2-x1), by default 1

    ``v``: -0.5 (default) or 1 < q < 0

    ``c1, c2``: filter constants (default c1: 8, c2: 2)

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

    if q > 0:
        raise ValueError(
            "This algorithm works right now only for semi integration, i.e. (q<0)"
        )

    if q != -0.5:
        print("\nWarning: algorithm is only tested for v=0.5 and v=-0.5.")
        print("         Other values for v are not verified!\n")

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
