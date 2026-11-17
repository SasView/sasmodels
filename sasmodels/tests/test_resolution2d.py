"""Analytic checks for oriented-sample slit resolution."""
import numpy as np
from numpy import pi, sqrt

from sasmodels.direct_model import DirectModel
from sasmodels.resolution2d import Slit2D


def test_slit2d_transverse_length():
    """Integrate a transverse Gaussian over the long slit, including zero width."""
    from scipy.special import erf

    q = np.array([0.001, 0.002, 0.005])
    length, decay = 0.1, 50.
    expected = sqrt(pi/decay)*erf(length*sqrt(decay))
    for slit_length in [length, np.full_like(q, length)]:
        for width in [0., 0.001, np.array([0.001, 0.002, 0.003])]:
            res = Slit2D(q, q_length=slit_length, q_width=width, accuracy='xhigh')
            _, qy = res.q_calc
            actual = res.apply(np.exp(-decay*qy**2))
            np.testing.assert_allclose(actual, expected, rtol=1e-5)


def test_slit2d_scan_width():
    """Match the analytic quadratic moment of the scan-direction Gaussian."""
    from scipy.special import ndtr

    q = np.array([0.01, 0.02, 0.03])
    width = np.array([0.001, 0.0015, 0.002])
    length = 0.1
    res = Slit2D(q, q_length=length, q_width=width,
                 q_calc=np.linspace(0.001, 0.04, 3901))
    qx, _ = res.q_calc
    actual = res.apply(qx**2)
    # The existing scan resolution truncates its Gaussian at -2.5 and +3 sigma.
    lower, upper = -2.5, 3.
    norm = ndtr(upper) - ndtr(lower)
    density_lower = np.exp(-lower**2/2)/sqrt(2*pi)
    density_upper = np.exp(-upper**2/2)/sqrt(2*pi)
    mean = (density_lower-density_upper)/norm
    second = 1 + (lower*density_lower-upper*density_upper)/norm
    expected = 2*length*(q**2 + 2*q*width*mean + width**2*second)
    np.testing.assert_allclose(actual, expected, rtol=1e-4)


def test_oriented_slit_interface():
    """Read slit dimensions from Data1D and integrate a real model analytically."""
    from scipy.special import erf

    from sasmodels.core import load_model
    from sasmodels.data import Data1D

    q = np.array([0.001, 0.002, 0.005])
    length, sigma = 0.1, 0.1
    data = Data1D(x=q)
    data.oriented = True
    data.dxl = np.full_like(q, length)
    data.dxw = np.zeros_like(q)
    calculator = DirectModel(data, load_model('gaussian_peak', dtype='double', platform='dll'))
    scale, background = 2., 0.03
    actual = calculator(peak_pos=0., sigma=sigma, scale=scale, background=background)
    integral = sqrt(2*pi)*sigma*erf(length/(sqrt(2)*sigma))
    expected = scale*np.exp(-q**2/(2*sigma**2))*integral + background
    # DirectModel uses the default low-resolution transverse integration grid.
    np.testing.assert_allclose(actual, expected, rtol=1e-3)
