#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
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

class GpuTriEllipse:
    PARS = {'scale':1, 'axisA':35, 'axisB':100, 'axisC':400, 'sldEll':1e-6, 'sldSolv':6.3e-6, 'background':0,
            'theta':0, 'phi':0, 'psi':0}

    PD_PARS = ['axisA', 'axisB', 'axisC', 'theta', 'phi', 'psi']

    def __init__(self, qx, qy, dtype='float32'):
        ctx,_queue = card()
        src, qx, qy = set_precision(open('Kernel/Kernel-TriaxialEllipse.cpp').read(), qx, qy, dtype=dtype)
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
        axisA, axisB, axisC, theta, phi, psi = \
            [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
             for base in GpuTriEllipse.PD_PARS]

        axisA.value, axisA.weight = axisA.get_weights(pars['axisA'], 0, 10000, True)
        axisB.value, axisB.weight = axisB.get_weights(pars['axisB'], 0, 10000, True)
        axisC.value, axisC.weight = axisC.get_weights(pars['axisC'], 0, 10000, True)
        theta.value, theta.weight = theta.get_weights(pars['theta'], -90, 180, False)
        phi.value, phi.weight = phi.get_weights(pars['phi'], -90, 180, False)
        psi.value, psi.weight = psi.get_weights(pars['psi'], -90, 180, False)

        sum, norm, norm_vol, vol = 0.0, 0.0, 0.0, 0.0
        size = len(theta.weight)
        sub = pars['sldEll'] - pars['sldSolv']

        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64
        for a in xrange(len(axisA.weight)):
            for b in xrange(len(axisB.weight)):
                for c in xrange(len(axisC.weight)):
                    for t in xrange(len(theta.weight)):
                        for i in xrange(len(phi.weight)):
                            for s in xrange(len(psi.weight)):
                                self.prg.TriaxialEllipseKernel(queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b,
                                            real(sub), real(pars['scale']), real(axisA.value[a]), real(axisB.value[b]),
                                            real(axisC.value[c]), real(phi.value[i]), real(theta.value[t]), real(psi.value[s]),
                                            real(axisA.weight[a]), real(axisB.weight[b]), real(axisC.weight[c]), real(psi.weight[s]),
                                            real(phi.weight[i]), real(theta.weight[t]), np.uint32(self.qx.size), np.uint32(size))
                                cl.enqueue_copy(queue, self.res, self.res_b)
                                sum += self.res

                                vol += axisA.weight[a]*axisB.weight[b]*axisC.weight[c]*axisA.value[a]*axisB.value[b]*axisC.value[c]
                                norm_vol += axisA.weight[a]*axisB.weight[b]*axisC.weight[c]
                                norm += axisA.weight[a]*axisB.weight[b]*axisC.weight[c]*theta.weight[t]*phi.weight[i]*psi.weight[s]

      #  if size > 1:
       #     norm /= asin(1.0)

        if vol != 0.0 and norm_vol != 0.0:
            sum *= norm_vol/vol

        return sum/norm + pars['background']










































