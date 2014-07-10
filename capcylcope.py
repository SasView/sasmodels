#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from math import asin
import pyopencl as cl
from weights import GaussianDispersion

class GpuCapCylinder(object):
    PARS = {'scale':1, 'rad_cyl':1, 'rad_cap':1, 'length':1, 'sld_capcyl':1e-6, 'sld_solv':0, 'background':0,
             'theta':0, 'phi':0}

    PD_PARS = ['rad_cyl', 'length', 'rad_cap', 'theta', 'phi']

    def __init__(self, qx, qy):

        self.qx = np.asarray(qx, np.float32)
        self.qy = np.asarray(qy, np.float32)
        #create context, queue, and build program
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)

        self.prg = cl.Program(self.ctx, open('Kernel-CapCyl.cpp').read()).build()

        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(self.ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)
        self.vol_i = float(0.0)
        self.vol_b = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.vol_i.nbytes)

    def eval(self, pars):

        rad_cyl,length,rad_cap,theta,phi = \
            [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
             for base in GpuCapCylinder.PD_PARS]

        rad_cyl.value, rad_cyl.weight = rad_cyl.get_weights(pars['rad_cyl'], 0, 1000, True)
        rad_cap.value, rad_cap.weight = rad_cap.get_weights(pars['rad_cap'], 0, 1000, True)
        length.value, length.weight = length.get_weights(pars['length'], 0, 1000, True)
        theta.value, theta.weight = theta.get_weights(pars['theta'], -90, 180, False)
        phi.value, phi.weight = phi.get_weights(pars['phi'], -90, 180, False)

        sum, norm, norm_vol, vol = 0.0, 0.0, 0.0, 0.0
        size = len(theta.weight)
        sub = pars['sld_capcyl']-np.float32(['sld_solv'])

        for i in xrange(len(rad_cyl.weight)):
            for m in xrange(len(rad_cap.weight)):
                for j in xrange(len(length.weight)):
                    for k in xrange(len(theta.weight)):
                        for l in xrange(len(phi.weight)):

                            self.prg.CapCylinderKernel(self.queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b,
                                        self.vol_b, np.float32(rad_cyl.value[i]), np.float32(rad_cap.value[m]), np.float32(length.value[j]),
                                        np.float32(theta.value[k]), np.float32(phi.value[l]), np.float32(sub), np.float32(pars['scale']),
                                        np.float32(phi.weight[l]), np.float32(theta.weight[k]), np.float32(rad_cap.weight[m]),
                                        np.float32(rad_cyl.weight[i]), np.float32(length.weight[j]), np.uint32(self.qx.size), np.uint32(size))

                            cl.enqueue_copy(self.queue, self.res, self.res_b)
                            cl.enqueue_copy(self.queue, self.vol_i, self.vol_b)

                            sum += self.res
                            vol += rad_cyl.weight[i]*length.weight[j]*rad_cap.weight[m]*self.vol_i
                            norm_vol += rad_cyl.weight[i]*length.weight[j]*rad_cap.weight[m]
                            norm += rad_cyl.weight[i]*length.weight[j]*rad_cap.weight[m]*theta.weight[k]*phi.weight[l]

        if size > 1:
            norm /= asin(1.0)

        if vol != 0.0 and norm_vol != 0.0:
            sum *= norm_vol/vol

        return sum/norm + pars['background']








































