#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import pyopencl as cl
from weights import GaussianDispersion

class GpuCoreShellCylinder(object):
    PARS = {'scale':1, 'radius':1, 'thickness':1, 'length':1, 'core_sld':1e-6, 'shell_sld':1e-6, 'solvent_sld':0,
            'background':0, 'axis_theta':0, 'axis_phi':0}
    PD_PARS = ['radius', 'length', 'thickness', 'axis_phi', 'axis_theta']

    def __init__(self, qx, qy):
        self.qx = np.asarray(qx, np.float32)
        self.qy = np.asarray(qy, np.float32)
        #create context, queue, and build program
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        self.prg = cl.Program(self.ctx, open('Kernel-CoreShellCylinder.cpp').read()).build()

        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(self.ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):

        radius, length, thickness, axis_phi, axis_theta = [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
                                     for base in GpuCoreShellCylinder.PD_PARS]

        radius.value, radius.weight = radius.get_weights(pars['radius'], 0, 1000, True)
        length.value, length.weight = length.get_weights(pars['length'], 0, 1000, True)
        thickness.value, thickness.weight = thickness.get_weights(pars['thickness'], 0, 1000, True)
        axis_phi.value, axis_phi.weight = axis_phi.get_weights(pars['axis_phi'], -90, 180, False)
        axis_theta.value, axis_theta.weight = axis_theta.get_weights(pars['axis_theta'], -90, 180, False)

        sum, norm, norm_vol, vol = 0.0, 0.0, 0.0, 0.0
        size = len(axis_theta.weight)

        for i in xrange(len(radius.weight)):
            for j in xrange(len(length.weight)):
                for k in xrange(len(axis_theta.weight)):
                    for l in xrange(len(axis_phi.weight)):
                        for f in xrange(len(thickness.weight)):

                            self.prg.CoreShellCylinderKernel(self.queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b,
                                    np.float32(axis_theta.value[k]), np.float32(axis_phi.value[l]), np.float32(thickness.value[f]),
                                    np.float32(length.value[j]), np.float32(radius.value[i]), np.float32(pars['scale']),
                                    np.float32(radius.weight[i]), np.float32(length.weight[j]), np.float32(thickness.weight[f]),
                                    np.float32(axis_theta.weight[k]), np.float32(axis_phi.weight[l]), np.float32(pars['core_sld']),
                                    np.float32(pars['shell_sld']), np.float32(pars['solvent_sld']),np.uint32(size),
                                    np.uint32(self.qx.size))
                            cl.enqueue_copy(self.queue, self.res, self.res_b)

                            sum += self.res
                            vol += radius.weight[i]*length.weight[j]*thickness.weight[f]*pow(radius.value[i]+thickness.value[f],2)\
                                   *(length.value[j]+2.0*thickness.value[f])
                            norm_vol += radius.weight[i]*length.weight[j]*thickness.weight[k]
                            norm += radius.weight[i]*length.weight[j]*thickness.weight[f]*axis_theta.weight[k]\
                                    *axis_phi.weight[l]

        if size>1:
            norm /= math.asin(1.0)
        if vol != 0.0 and norm_vol != 0.0:
            sum *= norm_vol/vol

        return sum/norm + pars['background']
