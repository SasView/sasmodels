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

class GpuCoreShellCylinder(object):
    PARS = {'scale':1, 'radius':1, 'thickness':1, 'length':1, 'core_sld':1e-6, 'shell_sld':-1e-6, 'solvent_sld':0,
            'background':0, 'axis_theta':0, 'axis_phi':0}
    PD_PARS = ['radius', 'length', 'thickness', 'axis_phi', 'axis_theta']

    def __init__(self, qx, qy, dtype='float32'):
        #create context, queue, and build program
        ctx,_queue = card()
        src, qx, qy = set_precision(open('NR_BessJ1.cpp').read()+"\n"+open('Kernel-CoreShellCylinder.cpp').read(), qx, qy, dtype=dtype)
        self.prg = cl.Program(ctx, src).build()
        self.qx, self.qy = qx, qy


        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(qx)

    def eval(self, pars):

        _ctx,queue = card()
        radius, length, thickness, axis_phi, axis_theta = [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
                                     for base in GpuCoreShellCylinder.PD_PARS]

        radius.value, radius.weight = radius.get_weights(pars['radius'], 0, 1000, True)
        length.value, length.weight = length.get_weights(pars['length'], 0, 1000, True)
        thickness.value, thickness.weight = thickness.get_weights(pars['thickness'], 0, 1000, True)
        axis_phi.value, axis_phi.weight = axis_phi.get_weights(pars['axis_phi'], -90, 180, False)
        axis_theta.value, axis_theta.weight = axis_theta.get_weights(pars['axis_theta'], -90, 180, False)

        sum, norm, norm_vol, vol = 0.0, 0.0, 0.0, 0.0
        size = len(axis_theta.weight)

        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64
        for r in xrange(len(radius.weight)):
            for l in xrange(len(length.weight)):
                for at in xrange(len(axis_theta.weight)):
                    for p in xrange(len(axis_phi.weight)):
                        for th in xrange(len(thickness.weight)):

                            self.prg.CoreShellCylinderKernel(queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b,
                                    real(axis_theta.value[at]), real(axis_phi.value[p]), real(thickness.value[th]),
                                    real(length.value[l]), real(radius.value[r]), real(pars['scale']),
                                    real(radius.weight[r]), real(length.weight[l]), real(thickness.weight[th]),
                                    real(axis_theta.weight[at]), real(axis_phi.weight[p]), real(pars['core_sld']),
                                    real(pars['shell_sld']), real(pars['solvent_sld']),np.uint32(size),
                                    np.uint32(self.qx.size))
                            cl.enqueue_copy(queue, self.res, self.res_b)

                            sum += self.res
                            vol += radius.weight[r]*length.weight[l]*thickness.weight[th]*pow(radius.value[r]+thickness.value[th],2)\
                                   *(length.value[l]+2.0*thickness.value[th])
                            norm_vol += radius.weight[r]*length.weight[l]*thickness.weight[th]
                            norm += radius.weight[r]*length.weight[l]*thickness.weight[th]*axis_theta.weight[at]\
                                    *axis_phi.weight[p]

        if size>1:
            norm /= math.asin(1.0)
        if vol != 0.0 and norm_vol != 0.0:
            sum *= norm_vol/vol

        return sum/norm + pars['background']
