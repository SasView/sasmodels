#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pyopencl as cl

from weights import GaussianDispersion
from Models.sasmodel import card
import hi


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
        src, qx, qy = set_precision(open('TEST-Kernel-Ellipse.cpp').read(), qx, qy, dtype=dtype)
        self.prg = cl.Program(ctx, src).build()
        self.qx, self.qy = qx, qy
        place = np.ascontiguousarray(hi.place, dtype=int)
        #buffers
        mf = cl.mem_flags
        self.place_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=place)
        self.qy_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.qx_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.res_b = cl.Buffer(ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):
    #b_n = radius_b # want, a_n = radius_a # want, etc
        ctx,queue = card()
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

        x = [radius_a.value, radius_a.weight, radius_b.value, radius_b.weight, axis_theta.value,
                 axis_theta.weight, axis_phi.value, axis_phi.weight]
        array = np.hstack(x)

        array_b = cl.Buffer(ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=array)

        self.prg.EllipsoidKernel(queue, self.qx.shape, None, self.qx_b, self.qy_b, self.place_b, array_b, self.res_b,
                                 real(pars['scale']), real(sub), np.uint32(self.qx.size), np.uint32(len(axis_theta.weight)))
        #copy result back from buffer
        cl.enqueue_copy(queue, self.res, self.res_b)


        a = open("answer.txt", "w")
        for x in xrange(len(self.res)):
            a.write(str(self.res))
            a.write("\n")


        return self.res+pars['background']





