#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import pyopencl as cl
from weights import GaussianDispersion


class GpuCylinder(object):
    PARS = {
        'scale':1,'radius':1,'length':1,'sldCyl':1e-6,'sldSolv':0,'background':0,
        'cyl_theta':0,'cyl_phi':0,
    }
    PD_PARS = ['radius', 'length', 'cyl_theta', 'cyl_phi']

    def __init__(self, qx, qy):

        self.qx = np.asarray(qx, np.float32)
        self.qy = np.asarray(qy, np.float32)
        #create context, queue, and build program
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        self.prg = cl.Program(self.ctx, open('Kernel-Cylinder.cpp').read()).build()

        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(self.ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):

        radius,length,cyl_theta,cyl_phi = \
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

        #Loop over radius, length, theta, phi weight points
        for i in xrange(len(radius.weight)):
            for j in xrange(len(length.weight)):
                for k in xrange(len(cyl_theta.weight)):
                    for l in xrange(len(cyl_phi.weight)):

                        self.prg.CylinderKernel(self.queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b, np.float32(sub),
                                           np.float32(radius.value[i]), np.float32(length.value[j]), np.float32(pars['scale']),
                                           np.float32(radius.weight[i]), np.float32(length.weight[j]), np.float32(cyl_theta.weight[k]),
                                           np.float32(cyl_phi.weight[l]), np.float32(cyl_theta.value[k]), np.float32(cyl_phi.value[l]),
                                           np.uint32(self.qx.size), np.uint32(size))
                        cl.enqueue_copy(self.queue, self.res, self.res_b)
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



























