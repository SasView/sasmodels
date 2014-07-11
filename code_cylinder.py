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

class GpuCylinder(object):
    PARS = {
        'scale':1,'radius':1,'length':1,'sldCyl':1e-6,'sldSolv':0,'background':0,
        'cyl_theta':0,'cyl_phi':0,
    }
    PD_PARS = ['radius', 'length', 'cyl_theta', 'cyl_phi']

    def __init__(self, qx, qy, dtype='float32'):

        #create context, queue, and build program
        ctx,_queue = card()
        src, qx, qy = set_precision(open('NR_BessJ1.cpp').read()+"\n"+open('Kernel-Cylinder.cpp').read(), qx, qy, dtype=dtype)
        self.prg = cl.Program(ctx, src).build()
        self.qx, self.qy = qx, qy

        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):

        _ctx,queue = card()
        radius, length, cyl_theta, cyl_phi = \
            [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
             for base in GpuCylinder.PD_PARS]

        #Get the weights for each
        radius.value, radius.weight = radius.get_weights(pars['radius'], 0, 1000, True)
        length.value, length.weight = length.get_weights(pars['length'], 0, 1000, True)
        cyl_theta.value, cyl_theta.weight = cyl_theta.get_weights(pars['cyl_theta'], -90, 180, False)
        cyl_phi.value, cyl_phi.weight = cyl_phi.get_weights(pars['cyl_phi'], -90, 180, False)

        #Perform the computation, with all weight points
        sum, norm, norm_vol, vol = 0.0, 0.0, 0.0, 0.0
        size = len(cyl_theta.weight)
        sub = pars['sldCyl'] - pars['sldSolv']

        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64
        #Loop over radius, length, theta, phi weight points
        for i in xrange(len(radius.weight)):
            for j in xrange(len(length.weight)):
                for k in xrange(len(cyl_theta.weight)):
                    for l in xrange(len(cyl_phi.weight)):


                        self.prg.CylinderKernel(queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b, real(sub),
                                           real(radius.value[i]), real(length.value[j]), real(pars['scale']),
                                           real(radius.weight[i]), real(length.weight[j]), real(cyl_theta.weight[k]),
                                           real(cyl_phi.weight[l]), real(cyl_theta.value[k]), real(cyl_phi.value[l]),
                                           np.uint32(self.qx.size), np.uint32(size))
                        cl.enqueue_copy(queue, self.res, self.res_b)
                        sum += self.res
                        vol += radius.weight[i]*length.weight[j]*pow(radius.value[i], 2)*length.value[j]
                        norm_vol += radius.weight[i]*length.weight[j]
                        norm += radius.weight[i]*length.weight[j]*cyl_theta.weight[k]*cyl_phi.weight[l]

        if size > 1:
            norm /= math.asin(1.0)
        if vol != 0.0 and norm_vol != 0.0:
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

    pars = CylinderParameters(scale=1, radius=64.1, length=266.96, sldCyl=.291e-6, sldSolv=5.77e-6, background=0,
                              cyl_theta=0, cyl_phi=0)
    t = time()
    result = GpuCylinder(qx, qy)
    result.x = result.cylinder_fit(pars, r_n=10, t_n=10, l_n=10, p_n=10, r_w=.1, t_w=.1, l_w=.1, p_w=.1, sigma=3)
    result.x = np.reshape(result.x, r_shape)
    tt = time()

    print("Time taken: %f" % (tt - t))

    plt.pcolormesh(result.x)
    plt.show()


if __name__=="__main__":
    demo()



























