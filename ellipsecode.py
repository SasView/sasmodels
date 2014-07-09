#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import pyopencl as cl
from weights import GaussianDispersion


class GpuEllipse(object):
    PARS = {
    'scale':1, 'radius_a':1, 'radius_b':1, 'sldEll':1e-6, 'sldSolv':0, 'background':0, 'axis_theta':0, 'axis_phi':0,
    }
    PD_PARS = ['radius_a', 'radius_b', 'axis_theta', 'axis_phi']
    def __init__(self, qx, qy):

        self.qx = np.asarray(qx, np.float32)
        self.qy = np.asarray(qy, np.float32)
        #create context, queue, and build program
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        self.prg = cl.Program(self.ctx, open('Kernel-Ellipse.cpp').read()).build()

        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(self.ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):
    #b_n = radius_b # want, a_n = radius_a # want, etc
        radius_a, radius_b, axis_theta, axis_phi = \
            [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
             for base in GpuEllipse.PD_PARS]

        radius_a.value, radius_a.weight = radius_a.get_weights(pars['radius_a'], 0, 1000, True)
        radius_b.value, radius_b.weight = radius_b.get_weights(pars['radius_b'], 0, 1000, True)
        axis_theta.value, axis_theta.weight = axis_theta.get_weights(pars['axis_theta'], -90, 180, False)
        axis_phi.value, axis_phi.weight = axis_phi.get_weights(pars['axis_phi'], -90, 180, False)

        #Perform the computation, with all weight points
        sum, norm, norm_vol, vol = 0.0, 0.0, 0.0, 0.0
        size = len(axis_theta.weight)
        sub =  pars['sldEll'] - pars['sldSolv']

        #Loop over radius weight points
        for i in xrange(len(radius_a.weight)):
            #Loop over length weight points
            for j in xrange(len(radius_b.weight)):
                #Average over theta distribution
                for k in xrange(len(axis_theta.weight)):
                    #Average over phi distribution
                    for l in xrange(len(axis_phi.weight)):
                        #call the kernel
                        self.prg.EllipsoidKernel(self.queue, self.qx.shape, None, np.float32(radius_a.weight[i]),
                                        np.float32(radius_b.weight[j]), np.float32(axis_theta.weight[k]),
                                        np.float32(axis_phi.weight[l]), np.float32(pars['scale']), np.float32(radius_a.value[i]),
                                        np.float32(radius_b.value[j]), np.float32(sub),np.float32(pars['background']),
                                        np.float32(axis_theta.value[k]), np.float32(axis_phi.value[l]), self.qx_b, self.qy_b,
                                        self.res_b, np.uint32(self.qx.size), np.uint32(len(axis_theta.weight)))
                        #copy result back from buffer
                        cl.enqueue_copy(self.queue, self.res, self.res_b)
                        sum += self.res
                        vol += radius_a.weight[i]*radius_b.weight[j]*pow(radius_b.value[j], 2)*radius_a.value[i]
                        norm_vol += radius_a.weight[i]*radius_b.weight[j]
                        norm += radius_a.weight[i]*radius_b.weight[j]*axis_theta.weight[k]*axis_phi.weight[l]
        # Averaging in theta needs an extra normalization
        # factor to account for the sin(theta) term in the
        # integration (see documentation).
        if size > 1:
            norm /= math.asin(1.0)
        if vol.any() != 0.0 and norm_vol.any() != 0.0:
            sum *= norm_vol/vol

        return sum/norm+pars['background']


def demo():
    from time import time
    import matplotlib.pyplot as plt

    #create qx and qy evenly spaces
    qx = np.linspace(-.02, .02, 128)
    qy = np.linspace(-.02, .02, 128)
    qx, qy = np.meshgrid(qx, qy)

    #saved shape of qx
    r_shape = qx.shape
    #reshape for calculation; resize as float32
    qx = qx.flatten()
    qy = qy.flatten()

    #int main
    pars = EllipsoidParameters(.027, 60, 180, .297e-6, 5.773e-06, 4.9, 0, 90)

    t = time()
    result = GpuEllipse(qx, qy)
    result.x = result.ellipsoid_fit(qx, qy, pars, b_n=35, t_n=35, a_n=1, p_n=1, sigma=3, b_w=.1, t_w=.1, a_w=.1, p_w=.1)
    result.x = np.reshape(result.x, r_shape)
    tt = time()
    print("Time taken: %f" % (tt - t))

    plt.pcolormesh(result.x)
    plt.show()


if __name__ == "__main__":
    demo()




