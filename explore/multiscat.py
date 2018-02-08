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
"""

from __future__ import print_function, division

import argparse
import time

import numpy as np
from scipy.special import gamma

from sasmodels import core
from sasmodels import compare
from sasmodels import resolution2d
from sasmodels.resolution import Resolution, bin_edges
from sasmodels.data import empty_data1D, empty_data2D, plot_data
from sasmodels.direct_model import call_kernel


class MultipleScattering(Resolution):
    def __init__(self, qmax, nq, probability, is2d, window=2, power=0):
        self.qmax = qmax
        self.nq = nq
        self.power = power
        self.probability = probability
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

    def apply(self, theory, background=0.0, outfile="", plot=False):
        #t0 = time.time()
        if self.is2d:
            Iq_calc = theory
        else:
            q_corners = self.q_calc[0]
            Iq_calc = np.interp(self._qxy, q_corners, theory)
        Iq_calc = Iq_calc.reshape(self.nq, self.nq)
        #plotxy(Iq_calc); import pylab; pylab.figure()
        if self.power > 0:
            Iqxy = scattering_power(Iq_calc, self.power)
            powers = []
        else:
            Iqxy, powers = multiple_scattering(
                Iq_calc, self.probability,
                coverage=0.99,
                return_powers=(outfile!="") or plot,
                )
        #print("multiple scattering calc time", time.time()-t0)

        #plotxy(Iqxy); import pylab; pylab.figure()
        if self.is2d:
            if outfile and powers:
                data = np.vstack([Ipower.flatten() for Ipower in powers]).T
                np.savetxt(outfile + "_powers.txt", data)
            if outfile:
                data = np.vstack(Iq_calc).T
                np.savetxt(outfile + ".txt", data)
            if plot:
                import pylab
                plotxy(Iq_calc)
                pylab.title("single scattering")
                pylab.figure()

            return Iqxy + background
        else:
            # circular average, no anti-aliasing
            Iq = np.histogram(self._qxy, bins=self._edges, weights=Iqxy)[0]/self._norm
            if powers:
                Iq_powers = [np.histogram(self._qxy, bins=self._edges, weights=Ipower)[0]/self._norm
                             for Ipower in powers]#  or powers[:5] for the first five

            if outfile:
                data = np.vstack([q_corners, theory, Iq]).T
                np.savetxt(outfile + ".txt", data)
            if outfile and powers:
                # circular average, no anti-aliasing for individual powers
                data = np.vstack([q_corners] + Iq_powers).T
                np.savetxt(outfile + "_powers.txt", data)
            if plot:
                import pylab
                plotxy(Iqxy)
                pylab.title("multiple scattering")
                pylab.figure()
            if plot and powers:
                import pylab
                L = -np.log(1-self.probability)
                pylab.loglog(q_corners, Iq, label="total %g"%self.probability)
                for n, Ipower in enumerate(Iq_powers):
                    k = n+1
                    w = L**(k-1)/gamma(k+1)
                    pylab.loglog(q_corners, w*Ipower, label="scattering**%d"%k)
                pylab.legend()
                pylab.figure()

            return q_corners, Iq + background

def scattering_power(Iq, n):
    """
    Calculate the nth scattering power as a distribution.  To get the
    weighted contribution, scale by $\lambda^k e^{-\lambda}/k!$.
    """
    scale = np.sum(Iq)
    F = _forward_fft(Iq/scale)
    result = _inverse_fft(F**n)
    return result

def multiple_scattering(Iq, p, coverage=0.99, return_powers=False):
    """
    Compute multiple scattering for I(q) given scattering probability p.

    Given a probability p of scattering with the thickness, the expected
    number of scattering events, $\lambda$ is $-\log(1 - p)$, giving a
    Poisson weighted sum of single, double, triple, etc. scattering patterns.
    The number of patterns used is based on coverage (default 99%).
    """
    L = -np.log(1-p)
    num_scatter = truncated_poisson_invcdf(coverage, L)

    # Compute multiple scattering via convolution.
    scale = np.sum(Iq)
    F = _forward_fft(Iq/scale)
    coeffs = [L**(k-1)/gamma(k+1) for k in range(1, num_scatter+1)]
    multiple_scattering = F * np.polyval(coeffs[::-1], F)
    result = scale * _inverse_fft(multiple_scattering)

    if return_powers:
        powers = [scale * _inverse_fft(F**k) for k in range(1, num_scatter+1)]
    else:
        powers = []

    return result, powers

def truncated_poisson_invcdf(p, L):
    """
    Return smallest k such that cdf(k; L) > p from the truncated Poisson
    probability excluding k=0
    """
    # pmf(k; L) = L**k * exp(-L) / (k! * (1-exp(-L))
    cdf = 0
    pmf = -np.exp(-L) / np.expm1(-L)
    k = 0
    while cdf < p:
        k += 1
        pmf *= L/k
        cdf += pmf
    return k

def _forward_fft(Iq):
    # Prepare padded array and forward transform
    nq = Iq.shape[0]
    half_nq = nq//2
    frame = np.zeros((2*nq, 2*nq))
    frame[:half_nq, :half_nq] = Iq[half_nq:, half_nq:]
    frame[-half_nq:, :half_nq] = Iq[:half_nq, half_nq:]
    frame[:half_nq, -half_nq:] = Iq[half_nq:, :half_nq]
    frame[-half_nq:, -half_nq:] = Iq[:half_nq, :half_nq]
    fourier_frame = np.fft.fft2(frame)
    return fourier_frame

def _inverse_fft(fourier_frame):
    # Invert the transform and recover the transformed data
    nq = fourier_frame.shape[0]//2
    half_nq = nq//2
    frame = np.fft.ifft2(fourier_frame).real
    Iq = np.empty((nq, nq))
    Iq[half_nq:, half_nq:] = frame[:half_nq, :half_nq]
    Iq[:half_nq, half_nq:] = frame[-half_nq:, :half_nq]
    Iq[half_nq:, :half_nq] = frame[:half_nq, -half_nq:]
    Iq[:half_nq, :half_nq] = frame[-half_nq:, -half_nq:]
    return Iq

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
    parser.add_argument('-k', '--power', type=int, default=0, help="show pattern for nth scattering")
    parser.add_argument('-p', '--probability', type=float, default=0.1, help="scattering probability")
    parser.add_argument('-n', '--nq', type=int, default=1024, help='number of mesh points')
    parser.add_argument('-q', '--qmax', type=float, default=0.5, help='max q')
    parser.add_argument('-w', '--window', type=float, default=2.0, help='q calc = q max * window')
    parser.add_argument('-2', '--2d', dest='is2d', action='store_true', help='oriented sample')
    parser.add_argument('-s', '--seed', default=-1, help='random pars with given seed')
    parser.add_argument('-r', '--random', action='store_true', help='random pars with random seed')
    parser.add_argument('-o', '--outfile', type=str, default="", help='random pars with random seed')
    parser.add_argument('model', type=str, help='sas model name such as cylinder')
    parser.add_argument('pars', type=str, nargs='*', help='model parameters such as radius=30')
    opts = parser.parse_args()
    assert opts.nq%2 == 0, "require even # points"

    model = core.load_model(opts.model)
    pars = parse_pars(model, opts)
    res = MultipleScattering(opts.qmax, opts.nq, opts.probability, opts.is2d,
                             window=opts.window, power=opts.power)
    kernel = model.make_kernel(res.q_calc)
    #print(pars)
    bg = pars.get('background', 0.0)
    pars['background'] = 0.0
    Iq_calc = call_kernel(kernel, pars)
    Iq = res.apply(Iq_calc, background=bg, outfile=opts.outfile, plot=True)
    plotxy(Iq)
    import pylab
    if opts.power > 0:
        pylab.title('scattering power %d'%opts.power)
    else:
        pylab.title('multiple scattering with fraction %g'%opts.probability)
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
