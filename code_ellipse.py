#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import pyopencl as cl
from weights import GaussianDispersion
from sasmodel import card

def set_precision(src, qx, qy, dtype):
    qx = np.ascontiguousarray(qx, dtype=dtype)
    qy = np.ascontiguousarray(qy, dtype=dtype)
    if np.dtype(dtype) == np.dtype('float32'):
        header = """\
#define real float
"""
    else:
        header = """\
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define real double
"""
    return header+src, qx, qy

class GpuEllipse(object):
    PARS = {
    'scale':1, 'radius_a':1, 'radius_b':1, 'sldEll':1e-6, 'sldSolv':0, 'background':0, 'axis_theta':0, 'axis_phi':0,
    }
    PD_PARS = ['radius_a', 'radius_b', 'axis_theta', 'axis_phi']

    def __init__(self, qx, qy, dtype='float32'):

        ctx,_queue = card()
        src, qx, qy = set_precision(open('Kernel-Ellipse.cpp').read(), qx, qy, dtype=dtype)
        self.prg = cl.Program(ctx, src).build()
        self.qx, self.qy = qx, qy

        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):
    #b_n = radius_b # want, a_n = radius_a # want, etc
        _ctx,queue = card()
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
        sub = pars['sldEll'] - pars['sldSolv']
        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64

        #Loop over radius weight points
        for i in xrange(len(radius_a.weight)):
            #Loop over length weight points
            for j in xrange(len(radius_b.weight)):
                #Average over theta distribution
                for k in xrange(len(axis_theta.weight)):
                    #Average over phi distribution
                    for l in xrange(len(axis_phi.weight)):
                        #call the kernel
                        self.prg.EllipsoidKernel(queue, self.qx.shape, None, real(radius_a.weight[i]),
                                        real(radius_b.weight[j]), real(axis_theta.weight[k]),
                                        real(axis_phi.weight[l]), real(pars['scale']), real(radius_a.value[i]),
                                        real(radius_b.value[j]), real(sub), real(axis_theta.value[k]),
                                        real(axis_phi.value[l]), self.qx_b, self.qy_b, self.res_b,
                                        np.uint32(self.qx.size), np.uint32(len(axis_theta.weight)))
                        #copy result back from buffer
                        cl.enqueue_copy(queue, self.res, self.res_b)
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




