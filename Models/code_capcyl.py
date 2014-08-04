#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from math import sqrt, fabs, atan
import pyopencl as cl

from weights import GaussianDispersion
from sasmodel import card, set_precision



class GpuCapCylinder(object):
    PARS = {'scale':1, 'rad_cyl':1, 'rad_cap':1, 'length':1, 'sld_capcyl':1e-6, 'sld_solv':0, 'background':0,
             'theta':0, 'phi':0}

    PD_PARS = ['rad_cyl', 'length', 'rad_cap', 'theta', 'phi']

    def __init__(self, qx, qy, dtype='float32'):

        #create context, queue, and build program
        ctx,_queue = card()
        trala = open('Kernel/NR_BessJ1.cpp').read()+"\n"+open('Kernel/Capcyl_Kfun.cpp').read()+"\n"+open('Kernel/Kernel-CapCyl.cpp').read()
        src, qx, qy = set_precision(trala, qx, qy, dtype=dtype)
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
        self.res[:] = 0
        cl.enqueue_copy(queue, self.res_b, self.res)
        rad_cyl,length,rad_cap,theta,phi = \
            [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
             for base in GpuCapCylinder.PD_PARS]

        rad_cyl.value, rad_cyl.weight = rad_cyl.get_weights(pars['rad_cyl'], 0, 10000, True)
        rad_cap.value, rad_cap.weight = rad_cap.get_weights(pars['rad_cap'], 0, 10000, True)
        length.value, length.weight = length.get_weights(pars['length'], 0, 10000, True)
        theta.value, theta.weight = theta.get_weights(pars['theta'], -90, 180, False)
        phi.value, phi.weight = phi.get_weights(pars['phi'], -90, 180, False)

        sum, norm, norm_vol, vol = 0.0, 0.0, 0.0, 0.0
        size = len(theta.weight)
        sub = pars['sld_capcyl']-pars['sld_solv']
        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64

        for i in xrange(len(rad_cyl.weight)):
            for m in xrange(len(rad_cap.weight)):
                for j in xrange(len(length.weight)):
                    for k in xrange(len(theta.weight)):
                        for l in xrange(len(phi.weight)):
                            hDist = -1.0*sqrt(fabs(rad_cap.value[m]*rad_cap.value[m]-rad_cyl.value[i]*rad_cyl.value[i]))
                            vol_i = 4.0*atan(1.0)*rad_cyl.value[i]*rad_cyl.value[i]*length.value[j]+2.0*4.0*atan(1.0)/3.0\
                                            *((rad_cap.value[m]-hDist)*(rad_cap.value[m]-hDist)*(2*rad_cap.value[m]+hDist))
                            self.prg.CapCylinderKernel(queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b,
                                        real(vol_i), real(hDist), real(rad_cyl.value[i]), real(rad_cap.value[m]), real(length.value[j]),
                                        real(theta.value[k]), real(phi.value[l]), real(sub), real(pars['scale']),
                                        real(phi.weight[l]), real(theta.weight[k]), real(rad_cap.weight[m]),
                                        real(rad_cyl.weight[i]), real(length.weight[j]), real(theta.weight[k]), np.uint32(self.qx.size), np.uint32(size))

                            vol += rad_cyl.weight[i]*length.weight[j]*rad_cap.weight[m]*vol_i
                            norm_vol += rad_cyl.weight[i]*length.weight[j]*rad_cap.weight[m]
                            norm += rad_cyl.weight[i]*length.weight[j]*rad_cap.weight[m]*theta.weight[k]*phi.weight[l]

        #if size > 1:
         #   norm /= asin(1.0)
        cl.enqueue_copy(queue, self.res, self.res_b)
        sum += self.res
        if vol != 0.0 and norm_vol != 0.0:
            sum *= norm_vol/vol

        return sum/norm + pars['background']








































