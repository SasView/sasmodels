"""
#This software was developed by the University of Tennessee as part of the
#Distributed Data Analysis of Neutron Scattering Experiments (DANSE)
#project funded by the US National Science Foundation. 
#See the license text in license.txt
"""
from __future__ import division

import numpy as np
from numpy import pi, cos, sin, sqrt

from .resolution import Resolution

## Singular point
SIGMA_ZERO = 1.0e-010
## Limit of how many sigmas to be covered for the Gaussian smearing
# default: 2.5 to cover 98.7% of Gaussian
NSIGMA = 3.0
## Defaults
NR = {'xhigh':10, 'high':5, 'med':5, 'low':3}
NPHI ={'xhigh':20, 'high':12, 'med':6, 'low':4}

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
        self.nsigma= nsigma
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
        self.index = index

        self.qx_data = data.qx_data[index]
        self.qy_data = data.qy_data[index]
        self.q_data = data.q_data[index]

        dqx = getattr(data, 'dqx_data', None)
        dqy = getattr(data, 'dqy_data', None)
        if dqx is not None and dqy is not None:
            # Here dqx and dqy mean dq_parr and dq_perp
            self.dqx_data = dqx[index]
            self.dqy_data = dqy[index]
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
            .reshape(nq, nbins).transpose().flatten()
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
            qx_res = ( (dqx*cos(dphi) + q_r) * cos(-q_phi) +
                           dqy*sin(dphi) * sin(-q_phi))
            qy_res = (-(dqx*cos(dphi) + q_r) * sin(-q_phi) +
                           dqy*sin(dphi) * cos(-q_phi))
        else:
            qx_res = qx +  dqx*cos(dphi)
            qy_res = qy +  dqy*sin(dphi)


        return qx_res, qy_res, weight_res

    def apply(self, theory):
        if self.q_calc_weights is not None:
            # TODO: interpolate rather than recomputing all the different qx,qy
            # Resolution needs to be applied
            nq, nbins = len(self.qx_data), self.nr * self.nphi
            ## Reshape into 2d array to use np weighted averaging
            theory = np.reshape(theory, (nbins, nq))
            ## Averaging with Gaussian weighting: normalization included.
            value =np.average(theory, axis=0, weights=self.q_calc_weights)
            ## Return the smeared values in the range of self.index
            return value
        else:
            return theory

"""
if __name__ == '__main__':
    ## Test w/ 2D linear function
    x = 0.001*np.arange(1, 11)
    dx = np.ones(len(x))*0.0003
    y = 0.001*np.arange(1, 11)
    dy = np.ones(len(x))*0.001
    z = np.ones(10)
    dz = sqrt(z)

    from sas.dataloader import Data2D
    #for i in range(10): print i, 0.001 + i*0.008/9.0
    #for i in range(100): print i, int(math.floor( (i/ (100/9.0)) ))
    out = Data2D()
    out.data = z
    out.qx_data = x
    out.qy_data = y
    out.dqx_data = dx
    out.dqy_data = dy
    out.q_data = sqrt(dx * dx + dy * dy)
    index = np.ones(len(x), dtype = bool)
    out.mask = index
    from sas.models.LineModel import LineModel
    model = LineModel()
    model.setParam("A", 0)

    smear = Smearer2D(out, model, index)
    #smear.set_accuracy('Xhigh')
    value = smear.get_value()
    ## All data are ones, so the smeared should also be ones.
    print "Data length =", len(value)
    print " 2D linear function, I = 0 + 1*qy"
    text = " Gaussian weighted averaging on a 2D linear function will "
    text += "provides the results same as without the averaging."
    print text
    print "qx_data", "qy_data", "I_nonsmear", "I_smeared"
    for ind in range(len(value)):
        print x[ind], y[ind], model.evalDistribution([x, y])[ind], value[ind]


if __name__ == '__main__':
    ## Another Test w/ constant function
    x = 0.001*np.arange(1,11)
    dx = np.ones(len(x))*0.001
    y = 0.001*np.arange(1,11)
    dy = np.ones(len(x))*0.001
    z = np.ones(10)
    dz = sqrt(z)

    from DataLoader import Data2D
    #for i in range(10): print i, 0.001 + i*0.008/9.0
    #for i in range(100): print i, int(math.floor( (i/ (100/9.0)) ))
    out = Data2D()
    out.data = z
    out.qx_data = x
    out.qy_data = y
    out.dqx_data = dx
    out.dqy_data = dy
    index = np.ones(len(x), dtype = bool)
    out.mask = index
    from sas.models.Constant import Constant
    model = Constant()

    value = Smearer2D(out,model,index).get_value()
    ## All data are ones, so the smeared values should also be ones.
    print "Data length =",len(value), ", Data=",value
"""
