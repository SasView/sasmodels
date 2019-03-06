#!/usr/bin/env python
r"""
Multiple scattering calculator

Calculate multiple scattering using 2D FFT convolution.

Usage:

    -p, --probability: the scattering probability
    -q, --qmax: that max q that you care about
    -w, --window: the extension window (q is calculated for qmax*window)
    -n, --nq: the number of mesh points (dq = qmax*window/nq)
    -r, --random: generate a random parameter set
    -2, --2d: perform the calculation for an oriented pattern
    model_name
    model_par=value ...

Assume the probability of scattering is $p$. After each scattering event,
$1-p$ neutrons will leave the system and go to the detector, and the remaining
$p$ will scatter again.

Let the scattering probability for $n$ scattering event at $q$ be $f_n(q)$,
where
.. math:: f_1(q) = \frac{I_1(q)}{\int I_1(q) dq}
for $I_1(q)$, the single scattering from the system. After two scattering
events, the scattering probability will be the convolution of the first
scattering and itself, or $f_2(q) = (f_1*f_1)(q)$.  After $n$ events it will be
$f_n(q) = (f_1 * \cdots * f_1)(q)$.  The total scattering is calculated
as the weighted sum of $f_k$, with weights following the Poisson distribution
.. math:: P(k; \lambda) = \frac{\lambda^k e^{-\lambda}}{k!}
for $\lambda$ determined by the total thickness divided by the mean
free path between scattering, giving
.. math::
    I(q) = \sum_{k=0}^\infty P(k; \lambda) f_k(q)
The $k=0$ term is ignored since it involves no scattering.
We cut the series when cumulative probability is less than cutoff $C=99\%$,
which is $\max n$ such that
.. math::
    \frac{\sum_{k=1}^n \frac{P(k; \lambda)}{1 - P(0; \lambda)} < C

Using the convolution theorem, where
$F = \mathcal{F}(f)$ is the Fourier transform,
.. math:: f * g = \mathcal{F}^{-1}\{\mathcal{F}\{f\} \cdot \mathcal{F}\{g\}\}
so
.. math:: f * \ldots * f = \mathcal{F}^{-1}\{ F^n \}
Since the Fourier transform is a linear operator, we can move the polynomial
expression for the convolution into the transform, giving
.. math::
    I(q) = \mathcal{F}^{-1}\left\{ \sum_{k=1}^{n} P(k; \lambda) F^k \right\}
In the dilute limit $L \rightarrow 0$ only the $k=1$ term is active,
and so
.. math::
    P(1; \lambda) = \lambda e{-\lambda} = \int I_1(q) dq
therefore we compute
.. math::
    I(q) = \int I_1(q) dq \mathcal{F}^{-1}\left\{
        \sum_{l=1}^{n} \frac{P(k; \lambda)}{P(1; \lambda))} F^k \right\}

For speed we may use the fast fourier transform with a power of two.
The resulting $I(q)$ will be linearly spaced and likely heavily oversampled.
The usual pinhole or slit resolution calculation can performed from these
calculated values.

By default single precision OpenCL is used for the calculation.  Set the
environment variable *SAS_OPENCL=none* to use double precision numpy FFT
instead.  The OpenCL versions is about 10x faster on an elderly Mac with
Intel HD 4000 graphics.  The single precision numerical artifacts don't
seem to seriously impact overall accuracy, though they look pretty bad.
"""

from __future__ import print_function, division

import argparse
import time

import numpy as np
from numpy import pi
from scipy.special import gamma

from sasmodels import core
from sasmodels import compare
from sasmodels.resolution import Resolution, bin_edges
from sasmodels.direct_model import call_kernel
import sasmodels.kernelcl

# TODO: select fast and accurate fft library
# clFFT: https://github.com/clMathLibraries/clFFT (AMD's OpenCL)
# - gpyfft: https://github.com/geggo/gpyfft (wraps clFFT)
# - arrayfire: https://github.com/arrayfire (wraps clFFT and much more; +cuda)
# pyFFT: https://github.com/fjarri-attic/pyfft (based on Apple's OpenCL; +cuda)
# - Reikna: https://github.com/fjarri/reikna (evolved from pyfft)
# genFFT: https://software.intel.com/en-us/articles/genFFT (Intel's OpenCL)
# VexCL: https://vexcl.readthedocs.io/en/latest/ (c++ library)
# TODO: switch to fftw when opencl is not available

try:
    import pyfft.cl
    import pyopencl.array as cl_array
    HAVE_OPENCL = sasmodels.kernelcl.use_opencl()
except ImportError:
    HAVE_OPENCL = False
PRECISION = np.dtype('f' if HAVE_OPENCL else 'd')  # 'f' or 'd'
USE_FAST = True  # OpenCL faster, less accurate math

class ICalculator:
    """
    Multiple scattering calculator
    """
    def fft(self, Iq):
        """
        Compute the forward FFT for an image, real -> complex.
        """
        raise NotImplementedError()

    def ifft(self, Iq):
        """
        Compute the inverse FFT for an image, complex -> complex.
        """
        raise NotImplementedError()

    def mulitple_scattering(self, Iq):
        r"""
        Compute multiple scattering for I(q) given scattering probability p.

        Given a probability p of scattering with the thickness, the expected
        number of scattering events, $\lambda$ is $-\log(1 - p)$, giving a
        Poisson weighted sum of single, double, triple, etc. scattering patterns.
        The number of patterns used is based on coverage (default 99%).
        """
        raise NotImplementedError()

class NumpyCalculator(ICalculator):
    """
    Multiple scattering calculator using numpy fft.
    """
    def __init__(self, dims=None, dtype=PRECISION):
        self.dtype = dtype
        self.complex_dtype = np.dtype('F') if dtype == np.dtype('f') else np.dtype('D')

    def fft(self, Iq):
        #t0 = time.time()
        Iq = np.asarray(Iq, self.dtype)
        result = np.fft.fft2(Iq)
        #print("fft time", time.time()-t0)
        return result

    def ifft(self, fourier_frame):
        #t0 = time.time()
        fourier_frame = np.asarray(fourier_frame, self.complex_dtype)
        result = np.fft.ifft2(fourier_frame)
        #print("ifft time", time.time()-t0)
        return result

    def multiple_scattering(self, Iq, p, coverage=0.99):
        #t0 = time.time()
        coeffs = scattering_coeffs(p, coverage)
        poly = np.asarray(coeffs[::-1], dtype=self.dtype)
        scale = np.sum(Iq)
        frame = _forward_shift(Iq/scale, dtype=self.dtype)
        fourier_frame = np.fft.fft2(frame)
        convolved = fourier_frame * np.polyval(poly, fourier_frame)
        frame = np.fft.ifft2(convolved)
        result = scale * _inverse_shift(frame.real, dtype=self.dtype)
        #print("numpy multiscat time", time.time()-t0)
        return result

# polyval1(c, x) computes (...((c0 x + c1) x + c2) x ... + cn) x
# where c is an array of length *degree* and x is an array of
# complex values (type double2) of length *n*. 2-D arrays can of
# course be treated as 1-D arrays of length *nx* X *ny*.
# When compiling with sasmodels.kernelcl.compile_model the double precision
# types are converted to single precision as needed.  See the code in
# sasmodels.generate._convert_type for details.
POLYVAL1_KERNEL = """
kernel void polyval1(
    const int degree,
    global const double *coeff,
    const int n,
    global double2 *array)
{
    int index = get_global_id(0);
    if (index < n) {
        const double2 x = array[index];
        double2 total = coeff[0];
        for (int k=1; k < degree; k++) {
            total = fma(total, x, coeff[k]);
        }
        array[index] = total * x;
    }
}
"""

class OpenclCalculator(ICalculator):
    """
    Multiple scattering calculator using OpenCL via pyfft.
    """
    polyval1f = None
    polyval1d = None
    def __init__(self, dims, dtype=PRECISION):
        env = sasmodels.kernelcl.environment()
        context = env.get_context(dtype)
        if dtype == np.dtype('f'):
            if OpenclCalculator.polyval1f is None:
                program = sasmodels.kernelcl.compile_model(
                    context, POLYVAL1_KERNEL, dtype, fast=USE_FAST)
                # Assume context is always the same for a given dtype
                OpenclCalculator.polyval1f = program.polyval1
            self.dtype = dtype
            self.complex_dtype = np.dtype('F')
            self.polyval1 = OpenclCalculator.polyval1f
        else:
            if OpenclCalculator.polyval1d is None:
                program = sasmodels.kernelcl.compile_model(
                    context, POLYVAL1_KERNEL, dtype, fast=False)
                # Assume context is always the same for a given dtype
                OpenclCalculator.polyval1d = program.polyval1
            self.dtype = dtype
            self.complex_type = np.dtype('D')
            self.polyval1 = OpenclCalculator.polyval1d
        self.queue = env.get_queue(dtype)
        self.plan = pyfft.cl.Plan(dims, queue=self.queue)

    def fft(self, Iq):
        # forward transform
        #t0 = time.time()
        data = np.asarray(Iq, self.complex_dtype)
        gpu_data = cl_array.to_device(self.queue, data)
        self.plan.execute(gpu_data.data)
        result = gpu_data.get()
        #print("fft time", time.time()-t0)
        return result

    def ifft(self, fourier_frame):
        # inverse transform
        #t0 = time.time()
        data = np.asarray(fourier_frame, self.complex_dtype)
        gpu_data = cl_array.to_device(self.queue, data)
        self.plan.execute(gpu_data.data, inverse=True)
        result = gpu_data.get()
        #print("ifft time", time.time()-t0)
        return result

    def multiple_scattering(self, Iq, p, coverage=0.99):
        #t0 = time.time()
        coeffs = scattering_coeffs(p, coverage)
        scale = np.sum(Iq)
        poly = np.asarray(coeffs[::-1], self.dtype)
        frame = _forward_shift(Iq/scale, dtype=self.complex_dtype)
        gpu_data = cl_array.to_device(self.queue, frame)
        gpu_poly = cl_array.to_device(self.queue, poly)
        self.plan.execute(gpu_data.data)
        degree, data_size= poly.shape[0], frame.shape[0]*frame.shape[1]
        self.polyval1(
            self.queue, [data_size], None,
            np.int32(degree), gpu_poly.data, np.int32(data_size), gpu_data.data)
        self.plan.execute(gpu_data.data, inverse=True)
        frame = gpu_data.get()
        #result = scale * _inverse_shift(frame.real, dtype=self.dtype)
        result = scale * _inverse_shift(frame.real, dtype=self.dtype)
        #print("OpenCL multiscat time", time.time()-t0)
        return result

Calculator = OpenclCalculator if HAVE_OPENCL else NumpyCalculator

def scattering_powers(Iq, n, dtype='f', transform=None):
    r"""
    Calculate the scattering powers up to n.

    This includes 1 even though it should just be Iq itself

    The frames are unweighted; to weight scale by $\lambda^k e^{-\lambda}/k!$.
    """
    if transform is None:
        n_x, n_y = Iq.shape
        transform = Calculator(dims=(n_x*2, n_y*2), dtype=dtype)
    scale = np.sum(Iq)
    frame = _forward_shift(Iq/scale, dtype=dtype)
    F = transform.fft(frame)
    powers = [scale * _inverse_shift(transform.ifft(F**(k+1)).real, dtype=dtype)
              for k in range(n)]
    return powers

def scattering_coeffs(p, coverage=0.99):
    r"""
    Return the coefficients of the scattering powers for transmission
    probability *p*.  This is just the corresponding values for the
    Poisson distribution for $\lambda = -\ln(1-p)$ such that
    $\sum_{k = 0 \ldots n} P(k; \lambda)$ is larger than *coverage*.
    """
    L = -np.log(1-p)
    num_scatter = truncated_poisson_invcdf(coverage, L)
    coeffs = [L**k/gamma(k+2) for k in range(num_scatter)]
    return coeffs

def truncated_poisson_invcdf(coverage, L):
    r"""
    Return smallest k such that cdf(k; L) > coverage from the truncated Poisson
    probability excluding k=0
    """
    # pmf(k; L) = L**k * exp(-L) / (k! * (1-exp(-L))
    cdf = 0
    pmf = -np.exp(-L) / np.expm1(-L)
    k = 0
    while cdf < coverage:
        k += 1
        pmf *= L/k
        cdf += pmf
    return k

def _forward_shift(Iq, dtype=PRECISION):
    # Prepare padded array and forward transform
    nq = Iq.shape[0]
    half_nq = nq//2
    frame = np.zeros((2*nq, 2*nq), dtype=dtype)
    frame[:half_nq, :half_nq] = Iq[half_nq:, half_nq:]
    frame[-half_nq:, :half_nq] = Iq[:half_nq, half_nq:]
    frame[:half_nq, -half_nq:] = Iq[half_nq:, :half_nq]
    frame[-half_nq:, -half_nq:] = Iq[:half_nq, :half_nq]
    return frame

def _inverse_shift(frame, dtype=PRECISION):
    # Invert the transform and recover the transformed data
    nq = frame.shape[0]//2
    half_nq = nq//2
    Iq = np.empty((nq, nq), dtype=dtype)
    Iq[half_nq:, half_nq:] = frame[:half_nq, :half_nq]
    Iq[:half_nq, half_nq:] = frame[-half_nq:, :half_nq]
    Iq[half_nq:, :half_nq] = frame[:half_nq, -half_nq:]
    Iq[:half_nq, :half_nq] = frame[-half_nq:, -half_nq:]
    return Iq


class MultipleScattering(Resolution):
    r"""
    Compute multiple scattering using Fourier convolution.

    The fourier steps are determined by *qmax*, the maximum $q$ value
    desired, *nq* the number of $q$ steps and *window*, the amount
    of padding around the circular convolution.  The $q$ spacing
    will be $\Delta q = 2 q_\mathrm{max} w / n_q$.  If *nq* is not
    given it will use $n_q = 2^k$ such that $\Delta q < q_\mathrm{min}$.

    *probability* is related to the expected number of scattering
    events in the sample $\lambda$ as $p = 1 - e^{-\lambda}$.
    *coverage* determines how many scattering steps to consider.  The
    default is 0.99, which sets $n$ such that $1 \ldots n$ covers 99%
    of the Poisson probability mass function.

    *is2d* is True then 2D scattering is used, otherwise it accepts
    and returns 1D scattering.

    *resolution* is the resolution function to apply after multiple
    scattering.  If present, then the resolution $q$ vectors will provide
    default values for *qmin*, *qmax* and *nq*.
    """
    def __init__(self, qmin=None, qmax=None, nq=None, window=2,
                 probability=None, coverage=0.99,
                 is2d=False, resolution=None,
                 dtype=PRECISION):
        # Infer qmin, qmax from instrument resolution calculator, if present
        if resolution is not None:
            is2d = hasattr(resolution, 'qx_data')
            if is2d:
                # 2D data
                if qmax is None:
                    qx_calc, qy_calc = resolution.q_calc
                    qmax = np.sqrt(np.max(qx_calc**2 + qy_calc**2))
                if qmin is None and nq is None:
                    qx, qy = resolution.data.x_bins, resolution.data.y_bins
                    if qx and qy:
                        dx = (np.max(qx) - np.min(qx)) / len(qx)
                        dy = (np.max(qy) - np.min(qy)) / len(qy)
                    else:
                        qx, qy = resolution.data.qx_data, resolution.data.qy_data
                        steps = np.sqrt(len(qx))
                        dx = (np.max(qx) - np.min(qx)) / steps
                        dy = (np.max(qy) - np.min(qy)) / steps
                    qmin = min(dx, dy)
            else:
                # 1D data
                if qmax is None:
                    qmax = np.max(resolution.q_calc)
                if qmin is None and nq is None:
                    qmin = np.min(np.abs(resolution.q_calc))

        # estimate nq from qmin, qmax if not given explicitly
        q_range = qmax * window
        if nq is None:
            nq = 2**np.ceil(np.log2(q_range/qmin))
        nq = int(nq)
        # Compute available qmin based on nq
        qmin = 2*q_range / nq
        #print(nq)

        # remember input parameters
        self.qmax = qmax
        self.qmin = qmin
        self.nq = nq
        self.probability = 0. if probability is None else probability
        self.coverage = coverage
        self.is2d = is2d
        self.window = window
        self.resolution = resolution

        # Determine the q values to calculate
        q = np.linspace(-q_range, q_range, nq)
        qx, qy = np.meshgrid(q, q)
        if is2d:
            q_calc = (qx.flatten(), qy.flatten())
        else:
            # For 1-D patterns, compute q from the center to the corners and
            # interpolate from there into the individual pixels.  Given that
            # nq represents the points in [-qmax*windows, qmax*window],
            # then using 2*sqrt(2)*nq/2 will oversample the points in q by
            # a factor of two relative to the pixels.
            q_range_to_corner = np.sqrt(2.) * q_range
            nq_to_corner = 10*int(np.ceil(np.sqrt(2.) * nq))
            q_to_corner = np.linspace(0, q_range_to_corner, nq_to_corner+1)[1:]
            q_calc = (q_to_corner,)
            # Remember the q radii of the calculated points
            self._radius = np.sqrt(qx**2 + qy**2)
            #self._q = q_to_corner
        self._q_steps = q
        self.q_calc = q_calc

        # TODO: use cleaner data representation than that from sasview
        # Resolution function forwards underlying q data (for plotting, etc?)
        if is2d:
            if resolution is not None:
                # forward resolution function info to multiscattering
                self.qx_data = resolution.qx_data
                self.qy_data = resolution.qy_data
            else:
                # no underlying resolution function, but make it look like there is
                self.qx_data, self.qy_data = q_calc
        else:
            # 1-D radial profile is determined by the q values we need to
            # compute, either for the calculated q values for the resolution
            # function (if any) or for the raw q values desired
            self._q = np.linspace(qmin, qmax, nq//(2*window))
            self._edges = bin_edges(self._q)
            self._norm, _ = np.histogram(self._radius, bins=self._edges)
            if resolution is not None:
                self.q = resolution.q
            else:
                # no underlying resolution function, but make it look like there is
                self.q = self._q

        # Prepare the multiple scattering calculator (either numpy or OpenCL)
        self.transform = Calculator((2*nq, 2*nq), dtype=dtype)

        # Iq and Iqxy will be set during apply
        self.Iq = None # type: np.ndarray
        self.Iqxy = None # type: np.ndarray

        # Label probability as a fittable parameter, and give its external name
        # Note that the external name must be a valid python identifier, since
        # is will be set as an experiment attribute.
        self.fittable = {'probability': 'scattering_probability'}

    def apply(self, theory):
        if self.is2d:
            Iq_calc = theory
        else:
            Iq_calc = np.interp(self._radius, self.q_calc[0], theory)
        Iq_calc = Iq_calc.reshape(self.nq, self.nq)

        # CRUFT: don't need probability as a function anymore
        probability = self.probability() if callable(self.probability) else self.probability
        coverage = self.coverage
        #t0 = time.time()
        Iqxy = self.transform.multiple_scattering(Iq_calc, probability, coverage)
        #print("multiple scattering calc time", time.time()-t0)

        # remember the intermediate result in case we want to see it later
        self.Iqxy = Iqxy

        if self.is2d:
            if self.resolution is not None:
                Iqxy = self.resolution.apply(Iqxy)
            return Iqxy
        else:
            # remember the intermediate result in case we want to see it later
            Iq = self.radial_profile(Iqxy)
            self.Iq = Iq
            if self.resolution is not None:
                q = self._q
                Iq_res = np.interp(np.abs(self.resolution.q_calc), q, self.Iq)
                _ = """
                k = 6
                print("q theory", q[:k])
                print("Iq theory", theory[:k])
                print("interp NaN?", np.any(np.isnan(Iq_calc)))
                print("convolved NaN?", np.any(np.isnan(Iqxy)))
                print("Iq intgrated", self.Iq[:k])
                print("radius", self._radius[self.nq/2,self.nq/2:self.nq/2+k])
                print("q edges", self._edges[:k+1])
                print("Iq norm", self._norm[:k])
                print("res q", self.resolution.q_calc[:k])
                print("Iq res", Iq_res[:k])
                #print(Iq)
                #print(Iq_res)
                """
                Iq = self.resolution.apply(Iq_res)
            return Iq

    def radial_profile(self, Iqxy):
        """
        Compute that radial profile for the given Iqxy grid.  The grid should
        be defined as for
        """
        # circular average, no anti-aliasing
        Iq = np.histogram(self._radius, bins=self._edges, weights=Iqxy)[0]/self._norm
        return Iq


def annular_average(qxy, Iqxy, qbins):
    """
    Compute annular average of points in *Iqxy* at *qbins*.  The $q_x$, $q_y$
    coordinates for *Iqxy* are given in *qxy*.
    """
    qxy, Iqxy = qxy.flatten(), Iqxy.flatten()
    index = np.argsort(qxy)
    qxy, Iqxy = qxy[index], Iqxy[index]
    print(qxy.shape, Iqxy.shape, index.shape, qbins.shape)
    #values = rebin(np.vstack((0., qxy)), Iqxy, qbins)
    integral = np.cumsum(Iqxy)
    Io = np.diff(np.interp(qbins, qxy, integral, left=0.))
    # normalize by area of annulus
    # TODO: doesn't properly account for box.
    # Need to chop off chords that are greater than the box edges
    # (left, right, top and bottom), then add back in the corners
    # chopped off by both. https://en.wikipedia.org/wiki/Circular_segment
    norms = np.diff(pi*qbins**2)
    return Io/norms

def rebin(x, I, xo):
    """
    Rebin from edges *x*, bins I into edges *xo*.

    *x* and *xo* should be monotonically increasing.

    If *x* has duplicate values, then all corresponding values at I(x) will
    be effectively summed into the same bin.  If *xo* has duplicate values,
    the first bin will contain the entire contents and subsequent bins will
    contain zeros.
    """
    integral = np.cumsum(I)
    Io = np.diff(np.interp(xo, x[1:], integral, left=0.))
    return Io


def parse_pars(model, opts):
    # type: (ModelInfo, argparse.Namespace) -> Dict[str, float]
    """
    Parse par=val arguments from the command line.
    """

    seed = np.random.randint(1000000) if opts.random and opts.seed < 0 else opts.seed
    compare_opts = {
        'info': (model.info, model.info),
        'use_demo': False,
        'seed': seed,
        'mono': True,
        'magnetic': False,
        'values': opts.pars,
        #'show_pars': True,
        'show_pars': False,
        'is2d': opts.is2d,
    }
    # Note: sascomp allows comparison on a pair of models, so ignore the second.
    pars, _ = compare.parse_pars(compare_opts)
    return pars


def main():
    parser = argparse.ArgumentParser(
        description="Compute multiple scattering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    parser.add_argument('-p', '--probability', type=float, default=0.1,
                        help="scattering probability")
    parser.add_argument('-n', '--nq', type=int, default=1024,
                        help='number of mesh points')
    parser.add_argument('-q', '--qmax', type=float, default=0.5,
                        help='max q')
    parser.add_argument('-w', '--window', type=float, default=2.0,
                        help='q calc = q max * window')
    parser.add_argument('-2', '--2d', dest='is2d', action='store_true',
                        help='oriented sample')
    parser.add_argument('-s', '--seed', default=-1,
                        help='random pars with given seed')
    parser.add_argument('-r', '--random', action='store_true',
                        help='random pars with random seed')
    parser.add_argument('-o', '--outfile', type=str, default="",
                        help='random pars with random seed')
    parser.add_argument('model', type=str,
                        help='sas model name such as cylinder')
    parser.add_argument('pars', type=str, nargs='*',
                        help='model parameters such as radius=30')
    opts = parser.parse_args()
    assert opts.nq%2 == 0, "require even # points"

    model = core.load_model(opts.model)
    pars = parse_pars(model, opts)
    res = MultipleScattering(qmax=opts.qmax, nq=opts.nq, window=opts.window,
                             probability=opts.probability, is2d=opts.is2d)
    kernel = model.make_kernel(res.q_calc)
    #print(pars)
    bg = pars.get('background', 0.0)
    pars['background'] = 0.0
    theory = call_kernel(kernel, pars)
    Iq = res.apply(theory) + bg
    plot_and_save_powers(res, theory, Iq, outfile=opts.outfile, background=bg)

def plot_and_save_powers(res, theory, result, plot=True, outfile="", background=0.):
    import pylab
    probability, coverage = res.probability, res.coverage
    weights = scattering_coeffs(probability, coverage)

    # cribbed from MultipleScattering.apply
    if res.is2d:
        Iq_calc = theory
    else:
        Iq_calc = np.interp(res._radius, res.q_calc[0], theory)
    Iq_calc = Iq_calc.reshape(res.nq, res.nq)

    # Compute the scattering powers for 1, 2, ... n scattering events
    powers = scattering_powers(Iq_calc, len(weights))

    #plotxy(Iqxy); import pylab; pylab.figure()
    if res.is2d:
        if outfile:
            data = np.vstack([Ipower.flatten() for Ipower in powers]).T
            np.savetxt(outfile + "_powers.txt", data)
            data = np.vstack(Iq_calc).T
            np.savetxt(outfile + ".txt", data)
        if plot:
            plotxy((res._q_steps, res._q_steps), Iq_calc)
            pylab.title("single scattering F")
            for k, v in enumerate(powers[1:]):
                pylab.figure()
                plotxy((res._q_steps, res._q_steps), v+background)
                pylab.title("multiple scattering F^%d" % (k+2))
            pylab.figure()
            plotxy((res._q_steps, res._q_steps), res.Iqxy+background)
            pylab.title("total scattering for p=%g" % probability)
            if res.resolution is not None:
                pylab.figure()
                plotxy((res._q_steps, res._q_steps), result)
                pylab.title("total scattering with resolution")
    else:
        q = res._q
        Iq_powers = [res.radial_profile(Iqxy) for Iqxy in powers]
        if outfile:
            data = np.vstack([q, theory, res.Iq]).T
            np.savetxt(outfile + ".txt", data)
            # circular average, no anti-aliasing for individual powers
            data = np.vstack([q] + Iq_powers).T
            np.savetxt(outfile + "_powers.txt", data)
        if plot:
            # Plot 2D pattern for total scattering
            plotxy((res._q_steps, res._q_steps), res.Iqxy+background)
            pylab.title("total scattering for p=%g" % probability)
            pylab.figure()

            # Plot 1D pattern for partial scattering
            pylab.loglog(q, res.Iq+background, label="total for p=%g"%probability)
            if res.resolution is not None:
                pylab.loglog(q, result, label="total with dQ")
            #new_annulus = annular_average(res._radius, res.Iqxy, res._edges)
            #pylab.loglog(q, new_annulus+background, label="new total for p=%g"%probability)
            for n, (w, Ipower) in enumerate(zip(weights, Iq_powers)):
                pylab.loglog(q, w*Ipower+background, label="scattering^%d"%(n+1))
            pylab.legend()
            pylab.title('total scattering for p=%g' % probability)
    pylab.show()

def plotxy(q, Iq):
    import pylab
    # q is a tuple of (q,) or (qx, qy)
    if len(q) == 1:
        pylab.loglog(q[0], Iq)
    else:
        data = Iq.copy()
        data[Iq <= 0] = np.min(Iq[Iq > 0])/2
        pylab.imshow(np.log10(data))

if __name__ == "__main__":
    main()
