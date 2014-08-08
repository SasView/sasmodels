#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pyopencl as cl

from weights import GaussianDispersion
from sasmodel import card, set_precision

class GpuTriEllipse:
    PARS = {'scale':1, 'semi_axisA':35, 'semi_axisB':100, 'semi_axisC':400, 'sldEll':1e-6, 'sldSolv':6.3e-6, 'background':0,
            'axis_theta':0, 'axis_phi':0, 'axis_psi':0}

    PD_PARS = ['semi_axisA', 'semi_axisB', 'semi_axisC', 'axis_theta', 'axis_phi', 'axis_psi']

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
        self.res[:] = 0
        cl.enqueue_copy(queue, self.res_b, self.res)
        semi_axisA, semi_axisB, semi_axisC, axis_theta, axis_phi, axis_psi = \
            [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
             for base in GpuTriEllipse.PD_PARS]

        semi_axisA.value, semi_axisA.weight = semi_axisA.get_weights(pars['semi_axisA'], 0, 10000, True)
        semi_axisB.value, semi_axisB.weight = semi_axisB.get_weights(pars['semi_axisB'], 0, 10000, True)
        semi_axisC.value, semi_axisC.weight = semi_axisC.get_weights(pars['semi_axisC'], 0, 10000, True)
        axis_theta.value, axis_theta.weight = axis_theta.get_weights(pars['axis_theta'], -90, 180, False)
        axis_phi.value, axis_phi.weight = axis_phi.get_weights(pars['axis_phi'], -90, 180, False)
        axis_psi.value, axis_psi.weight = axis_psi.get_weights(pars['axis_psi'], -90, 180, False)

        sum, norm, norm_vol, vol = 0.0, 0.0, 0.0, 0.0
        size = len(axis_theta.weight)
        sub = pars['sldEll'] - pars['sldSolv']

        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64
        for a in xrange(len(semi_axisA.weight)):
            for b in xrange(len(semi_axisB.weight)):
                for c in xrange(len(semi_axisC.weight)):

                    vol += semi_axisA.weight[a]*semi_axisB.weight[b]*semi_axisC.weight[c]*semi_axisA.value[a]*semi_axisB.value[b]*semi_axisC.value[c]
                    norm_vol += semi_axisA.weight[a]*semi_axisB.weight[b]*semi_axisC.weight[c]

                    for t in xrange(len(axis_theta.weight)):
                        for i in xrange(len(axis_phi.weight)):
                            for s in xrange(len(axis_psi.weight)):
                                self.prg.TriaxialEllipseKernel(queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b,
                                            real(sub), real(pars['scale']), real(semi_axisA.value[a]), real(semi_axisB.value[b]),
                                            real(semi_axisC.value[c]), real(axis_phi.value[i]), real(axis_theta.value[t]), real(axis_psi.value[s]),
                                            real(semi_axisA.weight[a]), real(semi_axisB.weight[b]), real(semi_axisC.weight[c]), real(axis_psi.weight[s]),
                                            real(axis_phi.weight[i]), real(axis_theta.weight[t]), np.uint32(self.qx.size), np.uint32(size))

                                norm += semi_axisA.weight[a]*semi_axisB.weight[b]*semi_axisC.weight[c]*axis_theta.weight[t]*axis_phi.weight[i]*axis_psi.weight[s]


      #  if size > 1:
       #     norm /= asin(1.0)
        cl.enqueue_copy(queue, self.res, self.res_b)
        sum = self.res
        if vol != 0.0 and norm_vol != 0.0:
            sum *= norm_vol/vol

        return sum/norm + pars['background']










































