#!/usr/bin/env python
r"""
Multiple scattering calculator

Calculate multiple scattering using 2D FFT convolution.

Usage:

    -f, --fraction: the scattering probability
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

Let the scattering probability for a single scattering event at $q$ be $f(q)$,
where
.. math:: f(q) = p \frac{I_1(q)}{\int I_1(q) dq}
for $I_1(q)$, the single scattering from the system. After two scattering
events, the scattering will be the convolution of the first scattering
and itself, or $(f*f)(q)$.  After $n$ events it will be
$(f* \cdots * f)(q)$.

The total scattering after multiple events will then be the number that didn't
scatter the first time $(=1-p)$ plus the number that scattered only once
$(=(1-p)f)$ plus the number that scattered only twice $(=(1-p)(f*f))$, etc.,
so
.. math::
    I(q) = (1-p)\sum_{k=0}^{\infty} f^{*k}(q)
Since convolution is linear, the contributions from the individual terms
will fall off as $p^k$, so we can cut off the series
at $n = \lceil \ln C / \ln p \rceil$ for some small fraction $C$ (default is
$C = 0.001$).

Using the convolution theorem, where
$F = \mathcal{F}(f)$ is the Fourier transform,
.. math::

    f * g = \mathcal{F}^{-1}\{\mathcal{F}\{f\} \cdot \mathcal{F}\{g\}\}

so
.. math::

    f * \ldots * f = \mathcal{F}^{-1}\{ F^n \}

Since the Fourier transform is a linear operator, we can move the polynomial
expression for the convolution into the transform, giving
.. math::

    I(q) = \mathcal{F}^{-1}\left\{ (1-p) \sum_{k=0}^{n} F^k \right\}

We drop the transmission term, $k=0$, and rescale the result by the
total scattering $\int I_1(q) dq$.

For speed we use the fast fourier transform for a power of two.  The resulting
$I(q)$ will be linearly spaced and likely heavily oversampled.  The usual
pinhole or slit resolution calculation can performed from these calculated
values.
"""

from __future__ import print_function

import argparse
import time

import numpy as np

from sasmodels import core
from sasmodels import compare
from sasmodels import resolution2d
from sasmodels.resolution import Resolution, bin_edges
from sasmodels.data import empty_data1D, empty_data2D, plot_data
from sasmodels.direct_model import call_kernel


class MultipleScattering(Resolution):
    def __init__(self, qmax, nq, fraction, is2d, window=2, power=0):
        self.qmax = qmax
        self.nq = nq
        self.power = power
        self.fraction = fraction
        self.is2d = is2d
        self.window = window

        q_range = qmax * window
        q = np.linspace(-q_range, q_range, nq)
        qx, qy = np.meshgrid(q, q)

        if is2d:
            q_calc = [qx.flatten(), qy.flatten()]
        else:
            q_range_corners = np.sqrt(2.) * q_range
            nq_corners = int(np.sqrt(2.) * nq/2)
            q_corners = np.linspace(0, q_range_corners, nq_corners+1)[1:]
            q_calc = [q_corners]
            self._qxy = np.sqrt(qx**2 + qy**2)
            self._edges = bin_edges(q_corners)
            self._norm = np.histogram(self._qxy, bins=self._edges)[0]

        self.q_calc = q_calc
        self.q_range = q_range

    def apply(self, theory):
        t0 = time.time()
        if self.is2d:
            Iq_calc = theory
        else:
            q_corners = self.q_calc[0]
            Iq_calc = np.interp(self._qxy, q_corners, theory)
        Iq_calc = Iq_calc.reshape(self.nq, self.nq)
        #plotxy(Iq_calc); import pylab; pylab.figure()
        if self.power > 0:
            Iqxy = autoconv_n(Iq_calc, self.power)
        else:
            Iqxy = multiple_scattering(Iq_calc, self.fraction, cutoff=0.001)
        print("multiple scattering calc time", time.time()-t0)
        #plotxy(Iqxy); import pylab; pylab.figure()
        if self.is2d:
            if 1:
                import pylab
                plotxy(Iq_calc)
                pylab.title("single scattering")
                pylab.figure()

            return Iqxy
        else:
            # circular average, no anti-aliasing
            Iq = np.histogram(self._qxy, bins=self._edges, weights=Iqxy)[0]/self._norm

            if 1:
                import pylab
                pylab.loglog(q_corners, theory, label="single scattering")
                if self.power > 0:
                    label = "scattering power %d"%self.power
                else:
                    label = "scattering fraction %d"%self.fraction
                pylab.loglog(q_corners, Iq, label=label)
                pylab.legend()
                pylab.figure()
                return Iqxy
            return q_corners, Iq

def multiple_scattering(Iq_calc, frac, cutoff=0.001):
    #plotxy(Iq_calc)
    num_scatter = int(np.ceil(np.log(cutoff)/np.log(frac)))

    # Prepare padded array for transform
    nq = Iq_calc.shape[0]
    half_nq = nq//2
    frame = np.zeros((2*nq, 2*nq))
    frame[:half_nq, :half_nq] = Iq_calc[half_nq:, half_nq:]
    frame[-half_nq:, :half_nq] = Iq_calc[:half_nq, half_nq:]
    frame[:half_nq, -half_nq:] = Iq_calc[half_nq:, :half_nq]
    frame[-half_nq:, -half_nq:] = Iq_calc[:half_nq, :half_nq]
    #plotxy(frame)

    # Compute multiple scattering via convolution.
    scale = np.sum(Iq_calc)
    fourier_frame = np.fft.fft2(frame/scale)
    #plotxy(abs(frame))
    # total = (1-a)f + (1-a)af^2 + (1-a)a^2f^3 + ...
    #       = (1-a)f[1 + af + (af)^2 + (af)^3 + ...]
    multiple_scattering = (
        (1-frac)*fourier_frame
        *np.polyval(np.ones(num_scatter), frac*fourier_frame))
    conv_frame = scale*np.fft.ifft2(multiple_scattering).real

    # Recover the transformed data
    #plotxy(conv_frame)
    Iq_conv = np.empty((nq, nq))
    Iq_conv[half_nq:, half_nq:] = conv_frame[:half_nq, :half_nq]
    Iq_conv[:half_nq, half_nq:] = conv_frame[-half_nq:, :half_nq]
    Iq_conv[half_nq:, :half_nq] = conv_frame[:half_nq, -half_nq:]
    Iq_conv[:half_nq, :half_nq] = conv_frame[-half_nq:, -half_nq:]
    #plotxy(Iq_conv)
    return Iq_conv

def multiple_scattering_cl(Iq_calc, frac, cutoff=0.001):
    raise NotImplementedError("no support for opencl calculations at this time")

    import pyopencl as cl
    import pyopencl.array as cla
    from gpyfft.fft import FFT
    context = cl.create_some_context()
    queue = cl.CommandQueue(context)

    #plotxy(Iq_calc)
    num_scatter = int(np.ceil(np.log(cutoff)/np.log(frac)))

    # Prepare padded array for transform
    nq = Iq_calc.shape[0]
    half_nq = nq//2
    frame = np.zeros((2*nq, 2*nq), dtype='float32')
    frame[:half_nq, :half_nq] = Iq_calc[half_nq:, half_nq:]
    frame[-half_nq:, :half_nq] = Iq_calc[:half_nq, half_nq:]
    frame[:half_nq, -half_nq:] = Iq_calc[half_nq:, :half_nq]
    frame[-half_nq:, -half_nq:] = Iq_calc[:half_nq, :half_nq]
    #plotxy(frame)

    # Compute multiple scattering via convolution (OpenCL operations)
    frame_gpu = cla.to_device(queue, frame)
    fourier_frame_gpu = cla.zeros(frame.shape, dtype='complex64')
    scale = frame_gpu.sum()
    frame_gpu /= scale
    transform = FFT(context, queue, frame_gpu, fourier_frame_gpu, axes=(0,1))
    event, = transform.enqueue()
    event.wait()
    fourier_frame_gpu *= frac
    multiple_scattering_gpu = fourier_frame_gpu.copy()
    for _ in range(num_scatter-1):
        multiple_scattering_gpu += 1
        multiple_scattering_gpu *= fourier_frame_gpu
    multiple_scattering_gpu *= (1 - frac)/frac
    transform = FFT(context, queue, multiple_scattering_gpu, frame_gpu, axes=(0,1))
    event, = transform.enqueue(forward=False)
    event.wait()
    conv_frame = cla.from_device(queue, frame_gpu)

    # Recover the transformed data
    #plotxy(conv_frame)
    Iq_conv = np.empty((nq, nq))
    Iq_conv[half_nq:, half_nq:] = conv_frame[:half_nq, :half_nq]
    Iq_conv[:half_nq, half_nq:] = conv_frame[-half_nq:, :half_nq]
    Iq_conv[half_nq:, :half_nq] = conv_frame[:half_nq, -half_nq:]
    Iq_conv[:half_nq, :half_nq] = conv_frame[-half_nq:, -half_nq:]
    #plotxy(Iq_conv)
    return Iq_conv

def autoconv_n(Iq_calc, power):
    # Compute multiple scattering via convolution.
    #plotxy(Iq_calc)
    scale = np.sum(Iq_calc)
    nq = Iq_calc.shape[0]
    frame = np.zeros((2*nq, 2*nq))
    half_nq = nq//2
    frame[:half_nq, :half_nq] = Iq_calc[half_nq:, half_nq:]
    frame[-half_nq:, :half_nq] = Iq_calc[:half_nq, half_nq:]
    frame[:half_nq, -half_nq:] = Iq_calc[half_nq:, :half_nq]
    frame[-half_nq:, -half_nq:] = Iq_calc[:half_nq, :half_nq]
    #plotxy(frame)
    fourier_frame = np.fft.fft2(frame/scale)
    #plotxy(abs(frame))
    fourier_frame = fourier_frame**power
    conv_frame = scale*np.fft.ifft2(fourier_frame).real
    #plotxy(conv_frame)
    Iq_conv = np.empty((nq, nq))
    Iq_conv[half_nq:, half_nq:] = conv_frame[:half_nq, :half_nq]
    Iq_conv[:half_nq, half_nq:] = conv_frame[-half_nq:, :half_nq]
    Iq_conv[half_nq:, :half_nq] = conv_frame[:half_nq, -half_nq:]
    Iq_conv[:half_nq, :half_nq] = conv_frame[-half_nq:, -half_nq:]
    #plotxy(Iq_conv)
    return Iq_conv

def parse_pars(model, opts):
    # type: (ModelInfo, argparse.Namespace) -> Dict[str, float]

    seed = np.random.randint(1000000) if opts.random and opts.seed < 0 else opts.seed
    compare_opts = {
        'info': (model.info, model.info),
        'use_demo': False,
        'seed': seed,
        'mono': True,
        'magnetic': False,
        'values': opts.pars,
        'show_pars': True,
        'is2d': opts.is2d,
    }
    pars, pars2 = compare.parse_pars(compare_opts)
    return pars


def main():
    parser = argparse.ArgumentParser(
        description="Compute multiple scattering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    parser.add_argument('-p', '--power', type=int, default=0, help="show pattern for nth scattering")
    parser.add_argument('-f', '--fraction', type=float, default=0.1, help="scattering fraction")
    parser.add_argument('-n', '--nq', type=int, default=1024, help='number of mesh points')
    parser.add_argument('-q', '--qmax', type=float, default=0.5, help='max q')
    parser.add_argument('-w', '--window', type=float, default=2.0, help='q calc = q max * window')
    parser.add_argument('-2', '--2d', dest='is2d', action='store_true', help='oriented sample')
    parser.add_argument('-s', '--seed', default=-1, help='random pars with given seed')
    parser.add_argument('-r', '--random', action='store_true', help='random pars with random seed')
    parser.add_argument('model', type=str, help='sas model name such as cylinder')
    parser.add_argument('pars', type=str, nargs='*', help='model parameters such as radius=30')
    opts = parser.parse_args()
    assert opts.nq%2 == 0, "require even # points"

    model = core.load_model(opts.model)
    pars = parse_pars(model, opts)
    res = MultipleScattering(opts.qmax, opts.nq, opts.fraction, opts.is2d,
                             window=opts.window, power=opts.power)
    kernel = model.make_kernel(res.q_calc)
    #print(pars)
    bg = pars.get('background', 0.0)
    pars['background'] = 0.0
    Iq_calc = call_kernel(kernel, pars)
    t0 = time.time()
    for i in range(10):
        Iq_calc = call_kernel(kernel, pars)
    print("single scattering calc time", (time.time()-t0)/10)
    Iq = res.apply(Iq_calc) + bg
    plotxy(Iq)
    import pylab
    if opts.power > 0:
        pylab.title('scattering power %d'%opts.power)
    else:
        pylab.title('multiple scattering with fraction %g'%opts.fraction)
    pylab.show()

def plotxy(Iq):
    import pylab
    if isinstance(Iq, tuple):
        q, Iq = Iq
        pylab.loglog(q, Iq)
    else:
        data = Iq+0.
        data[Iq <= 0] = np.min(Iq[Iq>0])/2
        pylab.imshow(np.log10(data))
    #import pylab; pylab.show()

if __name__ == "__main__":
    main()
