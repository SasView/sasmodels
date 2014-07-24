#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pyopencl as cl
from weights import GaussianDispersion

def set_precision(src, qx, qy, dtype):
    qx = np.ascontiguousarray(qx, dtype=dtype)
    qy = np.ascontiguousarray(qy, dtype=dtype)
    if dtype == 'double':
        header = """\
#define real float
"""
    else:
        header = """\
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define real double
"""
    return header+src, qx, qy


class GpuLamellar(object):
    PARS = {
        'scale':1, 'bi_thick':1, 'sld_bi':1e-6, 'sld_sol':0, 'background':0,
    }
    PD_PARS = {'bi_thick'}
    def __init__(self, qx, qy, dtype='float32'):

        #create context, queue, and build program
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        src,qx,qy = set_precision(open('Kernel-Lamellar.cpp').read(), qx, qy, dtype=dtype)
        self.prg = cl.Program(self.ctx, src).build()
        self.qx, self.qy = qx, qy

        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(self.ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):

        bi_thick = GaussianDispersion(int(pars['bi_thick_pd_n']), pars['bi_thick_pd'], pars['bi_thick_pd_nsigma'])
        bi_thick.value, bi_thick.weight = bi_thick.get_weights(pars['bi_thick'], 0, 10000, True)

        sum, norm = 0.0, 0.0
        sub = pars['sld_bi'] - pars['sld_sol']

        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64
        for i in xrange(len(bi_thick.weight)):
            self.prg.LamellarKernel(self.queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b, real(bi_thick.value[i]),
                                    real(pars['scale']), real(sub), np.uint32(self.qx.size))
            cl.enqueue_copy(self.queue, self.res, self.res_b)

            sum += bi_thick.weight[i]*self.res
            norm += bi_thick.weight[i]

        return sum/norm + pars['background']





























