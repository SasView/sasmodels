r"""
Special Functions
.................

This following standard C99 math functions are available:

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
        or inaccurate on some platforms. Use sas_erf, sas_erfc and sas_gamma
        instead (see below). Note: lgamma(x) has not yet been tested.

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

The following special functions and scattering calculations are defined.
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

    p1evl(x, c, n):
        Evaluate normalized polynomial $p(x) = x^n + \sum_{i=0}^{n-1} c_i x^i$
        using Horner's method so it is faster and more accurate.

        $c = \{c_{n-1}, c_{n-2} \ldots, c_0 \}$ is the table of coefficients,
        sorted from highest to lowest.

    sas_gamma(x):
        Gamma function $\text{sas_gamma}(x) = \Gamma(x)$.

        The standard math function, tgamma(x) is unstable for $x < 1$
        on some platforms.

    sas_gammaln(x):
        log gamma function sas_gammaln\ $(x) = \log \Gamma(|x|)$.

        The standard math function, lgamma(x), is incorrect for single
        precision on some platforms.

    sas_gammainc(a, x), sas_gammaincc(a, x):
        Incomplete gamma function
        sas_gammainc\ $(a, x) = \int_0^x t^{a-1}e^{-t}\,dt / \Gamma(a)$
        and complementary incomplete gamma function
        sas_gammaincc\ $(a, x) = \int_x^\infty t^{a-1}e^{-t}\,dt / \Gamma(a)$

    sas_erf(x), sas_erfc(x):
        Error function
        $\text{sas_erf}(x) = \frac{2}{\sqrt\pi}\int_0^x e^{-t^2}\,dt$
        and complementary error function
        $\text{sas_erfc}(x) = \frac{2}{\sqrt\pi}\int_x^{\infty} e^{-t^2}\,dt$.

        The standard math functions erf(x) and erfc(x) are slower and broken
        on some platforms.

    sas_J0(x):
        Bessel function of the first kind $\text{sas_J0}(x)=J_0(x)$ where
        $J_0(x) = \frac{1}{\pi}\int_0^\pi \cos(x\sin(\tau))\,d\tau$.

        The standard math function j0(x) is not available on all platforms.

    sas_J1(x):
        Bessel function of the first kind  $\text{sas_J1}(x)=J_1(x)$ where
        $J_1(x) = \frac{1}{\pi}\int_0^\pi \cos(\tau - x\sin(\tau))\,d\tau$.

        The standard math function j1(x) is not available on all platforms.

    sas_JN(n, x):
        Bessel function of the first kind and integer order $n$:
        $\text{sas_JN}(n, x)=J_n(x)$ where
        $J_n(x) = \frac{1}{\pi}\int_0^\pi \cos(n\tau - x\sin(\tau))\,d\tau$.
        If $n$ = 0 or 1, it uses sas_J0(x) or sas_J1(x), respectively.

        The standard math function jn(n, x) is not available on all platforms.

    sas_Si(x):
        Sine integral $\text{Si}(x) = \int_0^x \tfrac{\sin t}{t}\,dt$.

        This function uses Taylor series for small and large arguments:

        For large arguments,

        .. math::

            \text{Si}(x) \sim \frac{\pi}{2}
             - \frac{\cos(x)}{x}
            \left(1 - \frac{2!}{x^2} + \frac{4!}{x^4} - \frac{6!}{x^6} \right)
             - \frac{\sin(x)}{x}
            \left(\frac{1}{x} - \frac{3!}{x^3} + \frac{5!}{x^5} - \frac{7!}{x^7}\right)

        For small arguments,

        .. math::

           \text{Si}(x) \sim x
           - \frac{x^3}{3\times 3!} + \frac{x^5}{5 \times 5!} - \frac{x^7}{7 \times 7!}
           + \frac{x^9}{9\times 9!} - \frac{x^{11}}{11\times 11!}

    sas_3j1x_x(x):
        Spherical Bessel form
        $\text{sph_j1c}(x) = 3 j_1(x)/x = 3 (\sin(x) - x \cos(x))/x^3$,
        with a limiting value of 1 at $x=0$, where $j_1(x)$ is the spherical
        Bessel function of the first kind and first order.

        This function uses a Taylor series for small $x$ for numerical accuracy.


    sas_2J1x_x(x):
        Bessel form $\text{sas_J1c}(x) = 2 J_1(x)/x$, with a limiting value
        of 1 at $x=0$, where $J_1(x)$ is the Bessel function of first kind
        and first order.


    gauss76.n, gauss76.z[i], gauss76.w[i]:
        Points $z_i$ and weights $w_i$ for 76-point Gaussian quadrature, respectively,
        computing $\int_{-1}^1 f(z)\,dz \approx \sum_{i=1}^{76} w_i\,f(z_i)$.
        When translating the model to C, include 'lib/gauss76.c' in the source
        and use :code:`GAUSS_N`, :code:`GAUSS_Z`, and :code:`GAUSS_W`.

        Similar arrays are available in :code:`gauss20` for 20-point quadrature
        and :code:`gauss150.c` for 150-point quadrature. By using
        :code:`import gauss76 as gauss` it is easy to change the number of
        points in the integration.
"""
# pylint: disable=unused-import

import numpy as np

# Functions to add to our standard set
from numpy import degrees, radians

# C99 standard math library functions
from numpy import exp, log, power as pow, expm1, log1p, sqrt, cbrt
from numpy import sin, cos, tan, arcsin as asin, arccos as acos, arctan as atan
from numpy import sinh, cosh, tanh, arcsinh as asinh, arccosh as acosh, arctanh as atanh
from numpy import arctan2 as atan2
from numpy import fabs, fmin, fmax, trunc, rint
from numpy import pi, nan, inf
from scipy.special import gamma as sas_gamma
from scipy.special import gammaln as sas_gammaln
from scipy.special import gammainc as sas_gammainc
from scipy.special import gammaincc as sas_gammaincc
from scipy.special import erf as sas_erf
from scipy.special import erfc as sas_erfc
from scipy.special import j0 as sas_J0
from scipy.special import j1 as sas_J1
from scipy.special import jn as sas_JN

# erf, erfc, tgamma, lgamma  **do not use**

# C99 standard math constants
M_PI, M_PI_2, M_PI_4, M_SQRT1_2, M_E = np.pi, np.pi/2, np.pi/4, np.sqrt(0.5), np.e
NAN = nan
INFINITY = inf

# non-standard constants
M_PI_180, M_4PI_3 = M_PI/180, 4*M_PI/3

# can't do SINCOS in python; use "s, c = SINCOS(x)" instead
def SINCOS(x):
    """return sin(x), cos(x)"""
    return sin(x), cos(x)
sincos = SINCOS

def square(x):
    """return x^2"""
    return x*x

def cube(x):
    """return x^3"""
    return x*x*x

def sas_sinx_x(x):
    """return sin(x)/x"""
    from numpy import sinc as _sinc
    return _sinc(x/M_PI)

def powr(x, y):
    """return x^y for x>0"""
    return x**y
def pown(x, n):
    """return x^n for n integer"""
    return x**n

FLOAT_SIZE = 8

def polevl(x, c, n):
    """return p(x) for polynomial p of degree n-1 with coefficients c"""
    return np.polyval(c[:n], x)

def p1evl(x, c, n):
    """return x^n + p(x) for polynomial p of degree n-1 with coefficients c"""
    return np.polyval(np.hstack(([1.], c))[:n], x)

def sas_Si(x):
    """return Si(x)"""
    from scipy.special import sici
    return sici(x)[0]

def sas_j1(x):
    """return j1(x)"""
    if np.isscalar(x):
        retvalue = (sin(x) - x*cos(x))/x**2 if x != 0. else 0.
    else:
        with np.errstate(all='ignore'):
            retvalue = (sin(x) - x*cos(x))/x**2
        retvalue[x == 0.] = 0.
    return retvalue

def sas_3j1x_x(x):
    """return 3*j1(x)/x"""
    if np.isscalar(x):
        retvalue = 3*(sin(x) - x*cos(x))/x**3 if x != 0. else 1.
    else:
        with np.errstate(all='ignore'):
            retvalue = 3*(sin(x) - x*cos(x))/x**3
        retvalue[x == 0.] = 1.
    return retvalue

def sas_2J1x_x(x):
    """return 2*J1(x)/x"""
    if np.isscalar(x):
        retvalue = 2*sas_J1(x)/x if x != 0 else 1.
    else:
        with np.errstate(all='ignore'):
            retvalue = 2*sas_J1(x)/x
        retvalue[x == 0] = 1.
    return retvalue


# Gaussians
class Gauss:
    def __init__(self, w, z):
        self.n = len(w)
        self.w = w
        self.z = z

gauss20 = Gauss(
    w=np.array([
        .0176140071391521,
        .0406014298003869,
        .0626720483341091,
        .0832767415767047,
        .10193011981724,
        .118194531961518,
        .131688638449177,
        .142096109318382,
        .149172986472604,
        .152753387130726,
        .152753387130726,
        .149172986472604,
        .142096109318382,
        .131688638449177,
        .118194531961518,
        .10193011981724,
        .0832767415767047,
        .0626720483341091,
        .0406014298003869,
        .0176140071391521
    ]),
    z=np.array([
        -.993128599185095,
        -.963971927277914,
        -.912234428251326,
        -.839116971822219,
        -.746331906460151,
        -.636053680726515,
        -.510867001950827,
        -.37370608871542,
        -.227785851141645,
        -.076526521133497,
        .0765265211334973,
        .227785851141645,
        .37370608871542,
        .510867001950827,
        .636053680726515,
        .746331906460151,
        .839116971822219,
        .912234428251326,
        .963971927277914,
        .993128599185095
    ])
)

gauss76 = Gauss(
    w=np.array([
        .00126779163408536,		#0
        .00294910295364247,
        .00462793522803742,
        .00629918049732845,
        .00795984747723973,
        .00960710541471375,
        .0112381685696677,
        .0128502838475101,
        .0144407317482767,
        .0160068299122486,
        .0175459372914742,		#10
        .0190554584671906,
        .020532847967908,
        .0219756145344162,
        .0233813253070112,
        .0247476099206597,
        .026072164497986,
        .0273527555318275,
        .028587223650054,
        .029773487255905,
        .0309095460374916,		#20
        .0319934843404216,
        .0330234743977917,
        .0339977794120564,
        .0349147564835508,
        .0357728593807139,
        .0365706411473296,
        .0373067565423816,
        .0379799643084053,
        .0385891292645067,
        .0391332242205184,		#30
        .0396113317090621,
        .0400226455325968,
        .040366472122844,
        .0406422317102947,
        .0408494593018285,
        .040987805464794,
        .0410570369162294,
        .0410570369162294,
        .040987805464794,
        .0408494593018285,		#40
        .0406422317102947,
        .040366472122844,
        .0400226455325968,
        .0396113317090621,
        .0391332242205184,
        .0385891292645067,
        .0379799643084053,
        .0373067565423816,
        .0365706411473296,
        .0357728593807139,		#50
        .0349147564835508,
        .0339977794120564,
        .0330234743977917,
        .0319934843404216,
        .0309095460374916,
        .029773487255905,
        .028587223650054,
        .0273527555318275,
        .026072164497986,
        .0247476099206597,		#60
        .0233813253070112,
        .0219756145344162,
        .020532847967908,
        .0190554584671906,
        .0175459372914742,
        .0160068299122486,
        .0144407317482767,
        .0128502838475101,
        .0112381685696677,
        .00960710541471375,		#70
        .00795984747723973,
        .00629918049732845,
        .00462793522803742,
        .00294910295364247,
        .00126779163408536		#75 (indexed from 0)
    ]),
    z=np.array([
        -.999505948362153,		#0
        -.997397786355355,
        -.993608772723527,
        -.988144453359837,
        -.981013938975656,
        -.972229228520377,
        -.961805126758768,
        -.949759207710896,
        -.936111781934811,
        -.92088586125215,
        -.904107119545567,		#10
        -.885803849292083,
        -.866006913771982,
        -.844749694983342,
        -.822068037328975,
        -.7980001871612,
        -.77258672828181,
        -.74587051350361,
        -.717896592387704,
        -.688712135277641,
        -.658366353758143,		#20
        -.626910417672267,
        -.594397368836793,
        -.560882031601237,
        -.526420920401243,
        -.491072144462194,
        -.454895309813726,
        -.417951418780327,
        -.380302767117504,
        -.342012838966962,
        -.303146199807908,		#30
        -.263768387584994,
        -.223945802196474,
        -.183745593528914,
        -.143235548227268,
        -.102483975391227,
        -.0615595913906112,
        -.0205314039939986,
        .0205314039939986,
        .0615595913906112,
        .102483975391227,			#40
        .143235548227268,
        .183745593528914,
        .223945802196474,
        .263768387584994,
        .303146199807908,
        .342012838966962,
        .380302767117504,
        .417951418780327,
        .454895309813726,
        .491072144462194,		#50
        .526420920401243,
        .560882031601237,
        .594397368836793,
        .626910417672267,
        .658366353758143,
        .688712135277641,
        .717896592387704,
        .74587051350361,
        .77258672828181,
        .7980001871612,	#60
        .822068037328975,
        .844749694983342,
        .866006913771982,
        .885803849292083,
        .904107119545567,
        .92088586125215,
        .936111781934811,
        .949759207710896,
        .961805126758768,
        .972229228520377,		#70
        .981013938975656,
        .988144453359837,
        .993608772723527,
        .997397786355355,
        .999505948362153		#75
    ])
)

gauss150 = Gauss(
    z=np.array([
        -0.9998723404457334,
        -0.9993274305065947,
        -0.9983473449340834,
        -0.9969322929775997,
        -0.9950828645255290,
        -0.9927998590434373,
        -0.9900842691660192,
        -0.9869372772712794,
        -0.9833602541697529,
        -0.9793547582425894,
        -0.9749225346595943,
        -0.9700655145738374,
        -0.9647858142586956,
        -0.9590857341746905,
        -0.9529677579610971,
        -0.9464345513503147,
        -0.9394889610042837,
        -0.9321340132728527,
        -0.9243729128743136,
        -0.9162090414984952,
        -0.9076459563329236,
        -0.8986873885126239,
        -0.8893372414942055,
        -0.8795995893549102,
        -0.8694786750173527,
        -0.8589789084007133,
        -0.8481048644991847,
        -0.8368612813885015,
        -0.8252530581614230,
        -0.8132852527930605,
        -0.8009630799369827,
        -0.7882919086530552,
        -0.7752772600680049,
        -0.7619248049697269,
        -0.7482403613363824,
        -0.7342298918013638,
        -0.7198995010552305,
        -0.7052554331857488,
        -0.6903040689571928,
        -0.6750519230300931,
        -0.6595056411226444,
        -0.6436719971150083,
        -0.6275578900977726,
        -0.6111703413658551,
        -0.5945164913591590,
        -0.5776035965513142,
        -0.5604390262878617,
        -0.5430302595752546,
        -0.5253848818220803,
        -0.5075105815339176,
        -0.4894151469632753,
        -0.4711064627160663,
        -0.4525925063160997,
        -0.4338813447290861,
        -0.4149811308476706,
        -0.3959000999390257,
        -0.3766465660565522,
        -0.3572289184172501,
        -0.3376556177463400,
        -0.3179351925907259,
        -0.2980762356029071,
        -0.2780873997969574,
        -0.2579773947782034,
        -0.2377549829482451,
        -0.2174289756869712,
        -0.1970082295132342,
        -0.1765016422258567,
        -0.1559181490266516,
        -0.1352667186271445,
        -0.1145563493406956,
        -0.0937960651617229,
        -0.0729949118337358,
        -0.0521619529078925,
        -0.0313062657937972,
        -0.0104369378042598,
        0.0104369378042598,
        0.0313062657937972,
        0.0521619529078925,
        0.0729949118337358,
        0.0937960651617229,
        0.1145563493406956,
        0.1352667186271445,
        0.1559181490266516,
        0.1765016422258567,
        0.1970082295132342,
        0.2174289756869712,
        0.2377549829482451,
        0.2579773947782034,
        0.2780873997969574,
        0.2980762356029071,
        0.3179351925907259,
        0.3376556177463400,
        0.3572289184172501,
        0.3766465660565522,
        0.3959000999390257,
        0.4149811308476706,
        0.4338813447290861,
        0.4525925063160997,
        0.4711064627160663,
        0.4894151469632753,
        0.5075105815339176,
        0.5253848818220803,
        0.5430302595752546,
        0.5604390262878617,
        0.5776035965513142,
        0.5945164913591590,
        0.6111703413658551,
        0.6275578900977726,
        0.6436719971150083,
        0.6595056411226444,
        0.6750519230300931,
        0.6903040689571928,
        0.7052554331857488,
        0.7198995010552305,
        0.7342298918013638,
        0.7482403613363824,
        0.7619248049697269,
        0.7752772600680049,
        0.7882919086530552,
        0.8009630799369827,
        0.8132852527930605,
        0.8252530581614230,
        0.8368612813885015,
        0.8481048644991847,
        0.8589789084007133,
        0.8694786750173527,
        0.8795995893549102,
        0.8893372414942055,
        0.8986873885126239,
        0.9076459563329236,
        0.9162090414984952,
        0.9243729128743136,
        0.9321340132728527,
        0.9394889610042837,
        0.9464345513503147,
        0.9529677579610971,
        0.9590857341746905,
        0.9647858142586956,
        0.9700655145738374,
        0.9749225346595943,
        0.9793547582425894,
        0.9833602541697529,
        0.9869372772712794,
        0.9900842691660192,
        0.9927998590434373,
        0.9950828645255290,
        0.9969322929775997,
        0.9983473449340834,
        0.9993274305065947,
        0.9998723404457334
    ]),
    w=np.array([
        0.0003276086705538,
        0.0007624720924706,
        0.0011976474864367,
        0.0016323569986067,
        0.0020663664924131,
        0.0024994789888943,
        0.0029315036836558,
        0.0033622516236779,
        0.0037915348363451,
        0.0042191661429919,
        0.0046449591497966,
        0.0050687282939456,
        0.0054902889094487,
        0.0059094573005900,
        0.0063260508184704,
        0.0067398879387430,
        0.0071507883396855,
        0.0075585729801782,
        0.0079630641773633,
        0.0083640856838475,
        0.0087614627643580,
        0.0091550222717888,
        0.0095445927225849,
        0.0099300043714212,
        0.0103110892851360,
        0.0106876814158841,
        0.0110596166734735,
        0.0114267329968529,
        0.0117888704247183,
        0.0121458711652067,
        0.0124975796646449,
        0.0128438426753249,
        0.0131845093222756,
        0.0135194311690004,
        0.0138484622795371,
        0.0141714592928592,
        0.0144882814685445,
        0.0147987907597169,
        0.0151028518701744,
        0.0154003323133401,
        0.0156911024699895,
        0.0159750356447283,
        0.0162520081211971,
        0.0165218992159766,
        0.0167845913311726,
        0.0170399700056559,
        0.0172879239649355,
        0.0175283451696437,
        0.0177611288626114,
        0.0179861736145128,
        0.0182033813680609,
        0.0184126574807331,
        0.0186139107660094,
        0.0188070535331042,
        0.0189920016251754,
        0.0191686744559934,
        0.0193369950450545,
        0.0194968900511231,
        0.0196482898041878,
        0.0197911283358190,
        0.0199253434079123,
        0.0200508765398072,
        0.0201676730337687,
        0.0202756819988200,
        0.0203748563729175,
        0.0204651529434560,
        0.0205465323660984,
        0.0206189591819181,
        0.0206824018328499,
        0.0207368326754401,
        0.0207822279928917,
        0.0208185680053983,
        0.0208458368787627,
        0.0208640227312962,
        0.0208731176389954,
        0.0208731176389954,
        0.0208640227312962,
        0.0208458368787627,
        0.0208185680053983,
        0.0207822279928917,
        0.0207368326754401,
        0.0206824018328499,
        0.0206189591819181,
        0.0205465323660984,
        0.0204651529434560,
        0.0203748563729175,
        0.0202756819988200,
        0.0201676730337687,
        0.0200508765398072,
        0.0199253434079123,
        0.0197911283358190,
        0.0196482898041878,
        0.0194968900511231,
        0.0193369950450545,
        0.0191686744559934,
        0.0189920016251754,
        0.0188070535331042,
        0.0186139107660094,
        0.0184126574807331,
        0.0182033813680609,
        0.0179861736145128,
        0.0177611288626114,
        0.0175283451696437,
        0.0172879239649355,
        0.0170399700056559,
        0.0167845913311726,
        0.0165218992159766,
        0.0162520081211971,
        0.0159750356447283,
        0.0156911024699895,
        0.0154003323133401,
        0.0151028518701744,
        0.0147987907597169,
        0.0144882814685445,
        0.0141714592928592,
        0.0138484622795371,
        0.0135194311690004,
        0.0131845093222756,
        0.0128438426753249,
        0.0124975796646449,
        0.0121458711652067,
        0.0117888704247183,
        0.0114267329968529,
        0.0110596166734735,
        0.0106876814158841,
        0.0103110892851360,
        0.0099300043714212,
        0.0095445927225849,
        0.0091550222717888,
        0.0087614627643580,
        0.0083640856838475,
        0.0079630641773633,
        0.0075585729801782,
        0.0071507883396855,
        0.0067398879387430,
        0.0063260508184704,
        0.0059094573005900,
        0.0054902889094487,
        0.0050687282939456,
        0.0046449591497966,
        0.0042191661429919,
        0.0037915348363451,
        0.0033622516236779,
        0.0029315036836558,
        0.0024994789888943,
        0.0020663664924131,
        0.0016323569986067,
        0.0011976474864367,
        0.0007624720924706,
        0.0003276086705538
    ])
)
