.. _Special_Functions:

Special Functions
=================

The C code follows the C99 standard, with the usual math functions,
as defined in
`OpenCL <https://www.khronos.org/registry/cl/sdk/1.1/docs/man/xhtml/mathFunctions.html>`_.
These are automatically available for your models.

Standard functions and constants include:

    M_PI, M_PI_2, M_PI_4, M_SQRT1_2, M_E:
        $\pi$, $\pi/2$, $\pi/4$, $1/\sqrt{2}$ and Euler's constant $e$
    exp, log, pow(x,y), expm1, log1p, sqrt, cbrt:
        Power functions $e^x$, $\ln x$, $x^y$, $e^x - 1$, $\ln 1 + x$,
        $\sqrt{x}$, $\sqrt[3]{x}$. The functions expm1(x) and log1p(x)
        are accurate across all $x$, including $x$ very close to zero.
    sin, cos, tan, asin, acos, atan:
        Trigonometry functions and inverses, operating on radians.
    sinh, cosh, tanh, asinh, acosh, atanh:
        Hyperbolic trigonometry functions.
    atan2(y,x):
        Angle from the $x$\ -axis to the point $(x,y)$, which is equal to
        $\tan^{-1}(y/x)$ corrected for quadrant.  That is, if $x$ and $y$ are
        both negative, then atan2(y,x) returns a value in quadrant III where
        atan(y/x) would return a value in quadrant I. Similarly for
        quadrants II and IV when $x$ and $y$ have opposite sign.
    fabs(x), fmin(x,y), fmax(x,y), trunc, rint:
        Floating point functions.  rint(x) returns the nearest integer.
    NAN:
        NaN, Not a Number, $0/0$.  Use isnan(x) to test for NaN.  Note that
        you cannot use :code:`x == NAN` to test for NaN values since that
        will always return false.  NAN does not equal NAN!  The alternative,
        :code:`x != x` may fail if the compiler optimizes the test away.
    INFINITY:
        $\infty, 1/0$.  Use isinf(x) to test for infinity, or isfinite(x)
        to test for finite and not NaN.
    erf, erfc, tgamma, lgamma:  **do not use**
        Special functions that should be part of the standard, but are missing
        or inaccurate on some platforms. Use sas_erf, sas_erfc, sas_gamma
        and sas_lgamma instead (see below).

Some non-standard constants and functions are also provided:

    M_PI_180, M_4PI_3:
        $\frac{\pi}{180}$, $\frac{4\pi}{3}$
    SINCOS(x, s, c):
        Macro which sets s=sin(x) and c=cos(x). The variables *c* and *s*
        must be declared first.
    square(x):
        $x^2$
    cube(x):
        $x^3$
    clip(a, a_min, a_max):
        $\min(\max(a, a_\text{min}), a_\text{max})$, or NaN if $a$ is NaN.
    sas_sinx_x(x):
        $\sin(x)/x$, with limit $\sin(0)/0 = 1$.
    powr(x, y):
        $x^y$ for $x \ge 0$; this is faster than general $x^y$ on some GPUs.
    pown(x, n):
        $x^n$ for $n$ integer; this is faster than general $x^n$ on some GPUs.
    FLOAT_SIZE:
        The number of bytes in a floating point value.  Even though all
        variables are declared double, they may be converted to single
        precision float before running. If your algorithm depends on
        precision (which is not uncommon for numerical algorithms), use
        the following::

            #if FLOAT_SIZE>4
            ... code for double precision ...
            #else
            ... code for single precision ...
            #endif
    SAS_DOUBLE:
        A replacement for :code:`double` so that the declared variable will
        stay double precision; this should generally not be used since some
        graphics cards do not support double precision.  There is no provision
        for forcing a constant to stay double precision.

The following special functions and scattering calculations are defined in
`sasmodels/models/lib <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib>`_.
These functions have been tuned to be fast and numerically stable down
to $q=0$ even in single precision.  In some cases they work around bugs
which appear on some platforms but not others, so use them where needed.
Add the files listed in :code:`source = ["lib/file.c", ...]` to your *model.py*
file in the order given, otherwise these functions will not be available.

    polevl(x, c, n):
        Polynomial evaluation $p(x) = \sum_{i=0}^n c_i x^i$ using Horner's
        method so it is faster and more accurate.

        $c = \{c_n, c_{n-1}, \ldots, c_0 \}$ is the table of coefficients,
        sorted from highest to lowest.

        :code:`source = ["lib/polevl.c", ...]` (`link to code <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/polevl.c>`_)

    p1evl(x, c, n):
        Evaluation of normalized polynomial $p(x) = x^n + \sum_{i=0}^{n-1} c_i x^i$
        using Horner's method so it is faster and more accurate.

        $c = \{c_{n-1}, c_{n-2} \ldots, c_0 \}$ is the table of coefficients,
        sorted from highest to lowest.

        :code:`source = ["lib/polevl.c", ...]`
        (`polevl.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/polevl.c>`_)

    sas_gamma(x):
        Gamma function sas_gamma\ $(x) = \Gamma(x)$.

        The standard math function, tgamma(x), is unstable for $x < 1$
        on some platforms.

        :code:`source = ["lib/sas_gamma.c", ...]`
        (`sas_gamma.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_gamma.c>`_)

    sas_gammaln(x):
        log gamma function sas_gammaln\ $(x) = \log \Gamma(|x|)$.

        The standard math function, lgamma(x), is incorrect for single
        precision on some platforms.

        :code:`source = ["lib/sas_gammainc.c", ...]`
        (`sas_gammainc.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_gammainc.c>`_)

    sas_gammainc(a, x), sas_gammaincc(a, x):
        Incomplete gamma function
        sas_gammainc\ $(a, x) = \int_0^x t^{a-1}e^{-t}\,dt / \Gamma(a)$
        and complementary incomplete gamma function
        sas_gammaincc\ $(a, x) = \int_x^\infty t^{a-1}e^{-t}\,dt / \Gamma(a)$

        :code:`source = ["lib/sas_gammainc.c", ...]`
        (`sas_gammainc.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_gammainc.c>`_)

    sas_erf(x), sas_erfc(x):
        Error function
        sas_erf\ $(x) = \frac{2}{\sqrt\pi}\int_0^x e^{-t^2}\,dt$
        and complementary error function
        sas_erfc\ $(x) = \frac{2}{\sqrt\pi}\int_x^{\infty} e^{-t^2}\,dt$.

        The standard math functions erf(x) and erfc(x) are slower and broken
        on some platforms.

        :code:`source = ["lib/polevl.c", "lib/sas_erf.c", ...]`
        (`sas_erf.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_erf.c>`_)

    sas_J0(x):
        Bessel function of the first kind sas_J0\ $(x)=J_0(x)$ where
        $J_0(x) = \frac{1}{\pi}\int_0^\pi \cos(x\sin(\tau))\,d\tau$.

        The standard math function j0(x) is not available on all platforms.

        :code:`source = ["lib/polevl.c", "lib/sas_J0.c", ...]`
        (`sas_J0.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_J0.c>`_)

    sas_J1(x):
        Bessel function of the first kind  sas_J1\ $(x)=J_1(x)$ where
        $J_1(x) = \frac{1}{\pi}\int_0^\pi \cos(\tau - x\sin(\tau))\,d\tau$.

        The standard math function j1(x) is not available on all platforms.

        :code:`source = ["lib/polevl.c", "lib/sas_J1.c", ...]`
        (`sas_J1.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_J1.c>`_)

    sas_JN(n, x):
        Bessel function of the first kind and integer order $n$,
        sas_JN\ $(n, x) =J_n(x)$ where
        $J_n(x) = \frac{1}{\pi}\int_0^\pi \cos(n\tau - x\sin(\tau))\,d\tau$.
        If $n$ = 0 or 1, it uses sas_J0($x$) or sas_J1($x$), respectively.

        Warning: JN(n,x) can be very inaccurate (0.1%) for x not in [0.1, 100].

        The standard math function jn(n, x) is not available on all platforms.

        :code:`source = ["lib/polevl.c", "lib/sas_J0.c", "lib/sas_J1.c", "lib/sas_JN.c", ...]`
        (`sas_JN.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_JN.c>`_)

    sas_Si(x):
        Sine integral Si\ $(x) = \int_0^x \tfrac{\sin t}{t}\,dt$.

        Warning: Si(x) can be very inaccurate (0.1%) for x in [0.1, 100].

        This function uses Taylor series for small and large arguments:

        For large arguments use the following Taylor series,

        .. math::

             \text{Si}(x) \sim \frac{\pi}{2}
             - \frac{\cos(x)}{x}\left(1 - \frac{2!}{x^2} + \frac{4!}{x^4} - \frac{6!}{x^6} \right)
             - \frac{\sin(x)}{x}\left(\frac{1}{x} - \frac{3!}{x^3} + \frac{5!}{x^5} - \frac{7!}{x^7}\right)

        For small arguments,

        .. math::

           \text{Si}(x) \sim x
           - \frac{x^3}{3\times 3!} + \frac{x^5}{5 \times 5!} - \frac{x^7}{7 \times 7!}
           + \frac{x^9}{9\times 9!} - \frac{x^{11}}{11\times 11!}

        :code:`source = ["lib/Si.c", ...]`
        (`Si.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_Si.c>`_)

    sas_3j1x_x(x):
        Spherical Bessel form
        sph_j1c\ $(x) = 3 j_1(x)/x = 3 (\sin(x) - x \cos(x))/x^3$,
        with a limiting value of 1 at $x=0$, where $j_1(x)$ is the spherical
        Bessel function of the first kind and first order.

        This function uses a Taylor series for small $x$ for numerical accuracy.

        :code:`source = ["lib/sas_3j1x_x.c", ...]`
        (`sas_3j1x_x.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_3j1x_x.c>`_)


    sas_2J1x_x(x):
        Bessel form sas_J1c\ $(x) = 2 J_1(x)/x$, with a limiting value
        of 1 at $x=0$, where $J_1(x)$ is the Bessel function of first kind
        and first order.

        :code:`source = ["lib/polevl.c", "lib/sas_J1.c", ...]`
        (`sas_J1.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/sas_J1.c>`_)


    Gauss76Z[i], Gauss76Wt[i]:
        Points $z_i$ and weights $w_i$ for 76-point Gaussian quadrature, respectively,
        computing $\int_{-1}^1 f(z)\,dz \approx \sum_{i=1}^{76} w_i\,f(z_i)$.

        Similar arrays are available in :code:`gauss20.c` for 20-point
        quadrature and in :code:`gauss150.c` for 150-point quadrature.
        The macros :code:`GAUSS_N`, :code:`GAUSS_Z` and :code:`GAUSS_W` are
        defined so that you can change the order of the integration by
        selecting an different source without touching the C code.

        :code:`source = ["lib/gauss76.c", ...]`
        (`gauss76.c <https://github.com/SasView/sasmodels/tree/master/sasmodels/models/lib/gauss76.c>`_)

Complex numbers
...............

A small complex number library is available using :code:`source = ["lib/cl_complex.h", ...]`. 
Numbers are defined as type cdouble, which can be passed to and returned  from
functions. Operations are as follows:

    |    declare:     cdouble x
    |    define:      cplx(real, imag)
    |    x.real:      creal(x)
    |    x.imag:      cimag(x)
    |    1j:          I
    |    real + x:    radd(real, x)
    |    real - x:    rsub(real, x)
    |    real * x:    rmul(real, x)
    |    real / x:    rdiv(real, x)
    |    x + real:    radd(real, x)
    |    x - real:    radd(-real, x)
    |    x * real:    rmul(real, x)
    |    x / real:    rmul(1.0/real, x)
    |    x + y:       cadd(x, y)
    |    x - y:       csub(x, y)
    |    x * y:       cmul(x, y)
    |    x / y:       cdiv(x, y)
    |    -x:          cneg(x)
    |    abs(x):      cabs(x)
    |    angle(x):    carg(x)
    |
    |    special functions:
    |        csqrt, cexp, cpow, clog, clog10
    |        csin, ccos, ctan
    |        csinh, ccosh, ctanh

.. _Python_Functions:

Python Functions
................

To ease the transition between python and C models, the :ref:`Special_Functions`
available in C models are reproduced in Python with the same names. All of
them are available for import from :mod:`sasmodels.special`. For example, you
can include the spherical form factor :code:`sas_3j1x_x` using
:code:`from sasmodels.special import sas_3j1x_x` in your python model. The
python functions operate on both scalars and arrays.

Unlike C models where the constant :code:`FLOAT_SIZE` is 4 or 8 depending on
whether the model is single or double precision,
:code:`sasmodels.special.FLOAT_SIZE` is fixed at 8. Python models should check
:code:`input.dtype == np.float32` or :code:`input.dtype == np.float64` if
needed.

:code:`SAS_DOUBLE` is for C type definitions; it is not available in Python.

:code:`sas_j1 = (sin(x) - x cos(x))/x**2` is an additional function that
is not available for C.

For Gauss-Lobatto integration points use
:code:`from sasmodels.special import gauss76 as gauss`, then access them
with attributes :code:`gauss.n`, :code:`gauss.z` and :code:`gauss.w`. This is
equivalent to the macros :code:`GAUSS_N`, :code:`GAUSS_Z` and :code:`GAUSS_W`
in C models.

The complex number support functions aren't needed in Python, and haven't been
implemented.

WARNING: The python functions do not yet have the same accuracy as the
corresponding C functions.
