#!/usr/bin/env python
# -*- coding: utf-8 -*-


import ctypes
from ctypes import c_int, c_double, c_void_p
import numpy as np
import pyopencl as cl
from pyopencl import mem_flags as mf

from weights import GaussianDispersion
from sasmodel import card, set_precision, set_precision_1d, tic, toc

class GpuEllipse(object):
    PARS = {
    'scale':1, 'radius_a':1, 'radius_b':1, 'sldEll':1e-6, 'sldSolv':0, 'background':0, 'axis_theta':0, 'axis_phi':0,
    }
    PD_PARS = ['radius_a', 'radius_b', 'axis_theta', 'axis_phi']

    def __init__(self, qx, qy, dtype='float32', cutoff=1e-5):
        ctx,_queue = card()
        src, qx, qy = set_precision(open('Kernel/Kernel-Ellipse_f.cpp').read(), qx, qy, dtype=dtype)
        self.prg = cl.Program(ctx, src).build()
        self.qx, self.qy = qx, qy
        self.cutoff = cutoff

        self.qx_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(ctx, mf.READ_WRITE, self.qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):
        tic()

        ctx, queue = card()

        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64
        loops, loop_lengths = make_loops(pars, dtype=self.qx.dtype)
        loops_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=loops)
        loops_l = cl.LocalMemory(len(loops.data))

        self.prg.EllipsoidKernel(queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b,
                                loops_b, loops_l, real(self.cutoff),
                                real(pars['scale']), real(pars['background']),
                                real(pars['sldEll']-pars['sldSolv']),
                                *[np.uint32(pn) for pn in loop_lengths])

        cl.enqueue_copy(queue, self.res, self.res_b)
        print toc()*1000, self.qx.shape[0]

        return self.res

class CpuCylinder(GpuEllipse):
    def __init__(self, qx, qy, dtype='float32', cutoff=1e-5):
        self.qx, self.qy = [np.ascontiguousarray(v,'d') for v in qx,qy]
        self.cutoff = cutoff
        self.res = np.empty_like(self.qx)
        self.dll = ctypes.CDLL('Kernel/ellipsoid.so')
        self.fn = self.dll['EllipsoidKernel']
        self.fn.argtypes = [
            c_void_p, c_void_p, c_void_p, c_int,
            c_void_p, c_double, c_double, c_double,
            c_int, c_int, c_int, c_int
        ]
    def eval(self, pars):
        loops, loop_lengths = make_loops(pars, dtype=self.qx.dtype)
        self.fn(self.qx.ctypes.data, self.qy.ctypes.data, self.res.ctypes.data, len(self.qx),
            loops.ctypes.data, self.cutoff, pars['scale'],
            pars['sldEll']-pars['sldSolv'], *loop_lengths)

        return self.res

def make_loops(pars, dtype='double'):
    # 0.2 ms on sparkle to form the final loops
    radius_a, radius_b, theta, phi = \
       [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
         for base in GpuEllipse.PD_PARS]
    parts = [
        radius_a.get_weights(pars['radius_a'], 0, 10000, True),
        radius_b.get_weights(pars['radius_b'], 0, 10000, True),
        theta.get_weights(pars['axis_theta'], -np.inf, np.inf, False),
        phi.get_weights(pars['axis_phi'], -np.inf, np.inf, False),
        ]
    # Make sure that weights are normalized to peaks at 1 so that
    # the tolerance term can be used properly on truncated distributions
    loops = np.hstack((v,w/w.max()) for v,w in parts)
    #loops = np.hstack(parts)
    loops = np.ascontiguousarray(loops.T, dtype).flatten()
    return loops, [len(p[0]) for p in parts]