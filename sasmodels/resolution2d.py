"""
#This software was developed by the University of Tennessee as part of the
#Distributed Data Analysis of Neutron Scattering Experiments (DANSE)
#project funded by the US National Science Foundation.
#See the license text in license.txt
"""
from __future__ import division

import numpy as np  # type: ignore
from numpy import pi, cos, sin, sqrt  # type: ignore

from . import resolution
from .resolution import Resolution

## Singular point
SIGMA_ZERO = 1.0e-010
## Limit of how many sigmas to be covered for the Gaussian smearing
# default: 2.5 to cover 98.7% of Gaussian
NSIGMA = 3.0
## Defaults
NR = {'xhigh':10, 'high':5, 'med':5, 'low':3}
NPHI = {'xhigh':20, 'high':12, 'med':6, 'low':4}

## Defaults
N_SLIT_PERP = {'xhigh':1000, 'high':500, 'med':200, 'low':50}
N_SLIT_PERP_DOC = ", ".join("%s=%d"%(name, value)
                            for value, name in
                            sorted((2*v+1, k) for k, v in N_SLIT_PERP.items()))

class Pinhole2D(Resolution):
    """
    Gaussian Q smearing class for SAS 2d data
    """

    def __init__(self, data=None, index=None,
                 nsigma=NSIGMA, accuracy='Low', coords='polar'):
        """
        Assumption: equally spaced bins in dq_r, dq_phi space.

        :param data: 2d data used to set the smearing parameters
        :param index: 1d array with len(data) to define the range
         of the calculation: elements are given as True or False
        :param nr: number of bins in dq_r-axis
        :param nphi: number of bins in dq_phi-axis
        :param coord: coordinates [string], 'polar' or 'cartesian'
        """
        ## Accuracy: Higher stands for more sampling points in both directions
        ## of r and phi.
        ## number of bins in r axis for over-sampling
        self.nr = NR[accuracy.lower()]
        ## number of bins in phi axis for over-sampling
        self.nphi = NPHI[accuracy.lower()]
        ## maximum nsigmas
        self.nsigma = nsigma
        self.coords = coords
        self._init_data(data, index)

    def _init_data(self, data, index):
        """
        Get qx_data, qy_data, dqx_data,dqy_data,
        and calculate phi_data=arctan(qx_data/qy_data)
        """
        # TODO: maybe don't need to hold copy of qx,qy,dqx,dqy,data,index
        # just need q_calc and weights
        self.data = data
        self.index = index if index is not None else slice(None)

        self.qx_data = data.qx_data[self.index]
        self.qy_data = data.qy_data[self.index]
        self.q_data = data.q_data[self.index]

        dqx = getattr(data, 'dqx_data', None)
        dqy = getattr(data, 'dqy_data', None)
        if dqx is not None and dqy is not None:
            # Here dqx and dqy mean dq_parr and dq_perp
            self.dqx_data = dqx[self.index]
            self.dqy_data = dqy[self.index]
            ## Remove singular points if exists
            self.dqx_data[self.dqx_data < SIGMA_ZERO] = SIGMA_ZERO
            self.dqy_data[self.dqy_data < SIGMA_ZERO] = SIGMA_ZERO
            qx_calc, qy_calc, weights = self._calc_res()
            self.q_calc = [qx_calc, qy_calc]
            self.q_calc_weights = weights
        else:
            # No resolution information
            self.dqx_data = self.dqy_data = None
            self.q_calc = [self.qx_data, self.qy_data]
            self.q_calc_weights = None

        #self.phi_data = np.arctan(self.qx_data / self.qy_data)

    def _calc_res(self):
        """
        Over sampling of r_nbins times phi_nbins, calculate Gaussian weights,
        then find smeared intensity
        """
        nr, nphi = self.nr, self.nphi
        # Total number of bins = # of bins
        nbins = nr * nphi
        # Number of bins in the dqr direction (polar coordinate of dqx and dqy)
        bin_size = self.nsigma / nr
        # in dq_r-direction times # of bins in dq_phi-direction
        # data length in the range of self.index
        nq = len(self.qx_data)

        # Mean values of dqr at each bins
        # starting from the half of bin size
        r = bin_size / 2.0 + np.arange(nr) * bin_size
        # mean values of qphi at each bines
        phi = np.arange(nphi)
        dphi = phi * 2.0 * pi / nphi
        dphi = dphi.repeat(nr)

        ## Transform to polar coordinate,
        #  and set dphi at each data points ; 1d array
        dphi = dphi.repeat(nq)
        q_phi = self.qy_data / self.qx_data

        # Starting angle is different between polar
        #  and cartesian coordinates.
        #if self.coords != 'polar':
        #    dphi += np.arctan( q_phi * self.dqx_data/ \
        #                  self.dqy_data).repeat(nbins).reshape(nq,\
        #                                nbins).transpose().flatten()

        # The angle (phi) of the original q point
        q_phi = np.arctan(q_phi).repeat(nbins)\
            .reshape([nq, nbins]).transpose().flatten()
        ## Find Gaussian weight for each dq bins: The weight depends only
        #  on r-direction (The integration may not need)
        weight_res = (np.exp(-0.5 * (r - bin_size / 2.0)**2)  -
                      np.exp(-0.5 * (r + bin_size / 2.0)**2))
        # No needs of normalization here.
        #weight_res /= np.sum(weight_res)
        weight_res = weight_res.repeat(nphi).reshape(nr, nphi)
        weight_res = weight_res.transpose().flatten()

        ## Set dr for all dq bins for averaging
        dr = r.repeat(nphi).reshape(nr, nphi).transpose().flatten()
        ## Set dqr for all data points
        dqx = np.outer(dr, self.dqx_data).flatten()
        dqy = np.outer(dr, self.dqy_data).flatten()

        qx = self.qx_data.repeat(nbins)\
            .reshape(nq, nbins).transpose().flatten()
        qy = self.qy_data.repeat(nbins)\
            .reshape(nq, nbins).transpose().flatten()

        # The polar needs rotation by -q_phi
        if self.coords == 'polar':
            q_r = sqrt(qx**2 + qy**2)
            qx_res = ((dqx*cos(dphi) + q_r) * cos(-q_phi)
                      + dqy*sin(dphi) * sin(-q_phi))
            qy_res = (-(dqx*cos(dphi) + q_r) * sin(-q_phi)
                      + dqy*sin(dphi) * cos(-q_phi))
        else:
            qx_res = qx + dqx*cos(dphi)
            qy_res = qy + dqy*sin(dphi)


        return qx_res, qy_res, weight_res

    def apply(self, theory):
        if self.q_calc_weights is not None:
            # TODO: interpolate rather than recomputing all the different qx,qy
            # Resolution needs to be applied
            nq, nbins = len(self.qx_data), self.nr * self.nphi
            ## Reshape into 2d array to use np weighted averaging
            theory = np.reshape(theory, (nbins, nq))
            ## Averaging with Gaussian weighting: normalization included.
            value = np.average(theory, axis=0, weights=self.q_calc_weights)
            ## Return the smeared values in the range of self.index
            return value
        else:
            return theory


class Slit2D(Resolution):
    """
    Slit aperture with resolution function on an oriented sample.

    *q* points at which the data is measured.

    *qx_width* slit width in qx

    *qy_width* slit height in qy; current implementation requires a fixed
    qy_width for all q points.

    *q_calc* is the list of q points to calculate, or None if this
    should be estimated from the *q* and *qx_width*.

    *accuracy* determines the number of *qy* points to compute for each *q*.
    The values are stored in sasmodels.resolution2d.N_SLIT_PERP.  The default
    values are: %s
    """
    __doc__ = __doc__%N_SLIT_PERP_DOC
    def __init__(self, q, qx_width, qy_width=0., q_calc=None, accuracy='low'):
        # Remember what q and width was used even though we won't need them
        # after the weight matrix is constructed
        self.q, self.qx_width, self.qy_width = q, qx_width, qy_width

        # Allow independent resolution on each qx point even though it is not
        # needed in practice.  Set qy_width to the maximum qy width.
        if np.isscalar(qx_width):
            qx_width = np.ones(len(q))*qx_width
        else:
            qx_width = np.asarray(qx_width)
        if not np.isscalar(qy_width):
            qy_width = np.max(qy_width)

        # Build grid of qx, qy points
        if q_calc is not None:
            qx_calc = np.sort(q_calc)
        else:
            qx_calc = resolution.pinhole_extend_q(q, qx_width, nsigma=3)
        qy_min, qy_max = np.log10(np.min(q)), np.log10(qy_width)
        qy_calc = np.logspace(qy_min, qy_max, N_SLIT_PERP[accuracy])
        qy_calc = np.hstack((-qy_calc[::-1], 0, qy_calc))
        self.q_calc = [v.flatten() for v in np.meshgrid(qx_calc, qy_calc)]
        self.qx_calc, self.qy_calc = qx_calc, qy_calc
        self.nx, self.ny = len(qx_calc), len(qy_calc)
        self.dy = 2*qy_width/self.ny

        # Build weight matrix for resolution integration
        if np.any(qx_width > 0):
            self.weights = resolution.pinhole_resolution(
                qx_calc, q, np.maximum(qx_width, resolution.MINIMUM_RESOLUTION))
        elif len(qx_calc) == len(q) and np.all(qx_calc == q):
            self.weights = None
        else:
            raise ValueError("Slit2D fails with q_calc != q")

    def apply(self, theory):
        Iq = np.trapz(theory.reshape(self.ny, self.nx), axis=0, x=self.qy_calc)
        if self.weights is not None:
            Iq = resolution.apply_resolution_matrix(self.weights, Iq)
        return Iq
