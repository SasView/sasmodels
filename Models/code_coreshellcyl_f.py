#!/usr/bin/env python
# -*- coding: utf-8 -*-


import ctypes
from ctypes import c_int, c_double, c_void_p
import numpy as np
import pyopencl as cl
from pyopencl import mem_flags as mf

from weights import GaussianDispersion
from sasmodel import card, set_precision, set_precision_1d, tic, toc

class GpuCoreShellCylinder(object):
    PARS = {'scale':1, 'radius':1, 'thickness':1, 'length':1, 'core_sld':1e-6, 'shell_sld':-1e-6, 'solvent_sld':0,
            'background':0, 'axis_theta':0, 'axis_phi':0}
    PD_PARS = ['radius', 'length', 'thickness', 'axis_phi', 'axis_theta']

    def __init__(self, qx, qy, dtype='float32', cutoff=1e-5):
        #create context, queue, and build program
        ctx,_queue = card()
        src, qx, qy = set_precision(open('Kernel/NR_BessJ1.cpp').read()+"\n"+open('Kernel/Kernel-CoreShellCylinder_f.cpp').read(), qx, qy, dtype=dtype)
        self.prg = cl.Program(ctx, src).build()
        self.qx, self.qy = qx, qy
        self.cutoff = cutoff

        self.qx_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(ctx, mf.READ_WRITE, self.qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):
        tic()

        ctx,queue = card()
        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64
        loops, loop_lengths = make_loops(pars, dtype=self.qx.dtype)
        loops_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=loops)
        loops_l = cl.LocalMemory(len(loops.data))

        self.prg.CoreShellCylinderKernel(queue, self.qx.shape, None,
                        self.qx_b, self.qy_b, self.res_b,
                        loops_b, loops_l, real(self.cutoff),
                        real(pars['scale']), real(pars['background']),
                        real(pars['shell_sld']-pars['solvent_sld']),
                        real(pars['core_sld']-pars['shell_sld']),
                        *[np.uint32(pn) for pn in loop_lengths])

        cl.enqueue_copy(queue, self.res, self.res_b)
        print toc()*1000, self.qx.shape[0]

        return self.res

class CPUCoreShellCylinder(GpuCoreShellCylinder):
    def __init__(self, qx, qy, dtype='float32', cutoff=1e-5):
        self.qx, self.qy = [np.ascontiguousarray(v,'d') for v in qx,qy]
        self.cutoff = cutoff
        self.res = np.empty_like(self.qx)
        self.dll = ctypes.CDLL('Kernel/coreshellcylinder.so')
        self.fn = self.dll['CoreShellCylinderKernel']
        self.fn.argtypes = [
            c_void_p, c_void_p, c_void_p, c_int,
            c_void_p, c_double, c_double, c_double, c_double, c_double,
            c_int, c_int, c_int, c_int, c_int
        ]
    def eval(self, pars):
        loops, loop_lengths = make_loops(pars, dtype=self.qx.dtype)
        self.fn(self.qx.ctypes.data, self.qy.ctypes.data, self.res.ctypes.data, len(self.qx),
            loops.ctypes.data, self.cutoff, pars['scale'], pars['background'],
            pars['shell_sld']-pars['solvent_sld'], pars['core_sld']-pars['shell-sld'], *loop_lengths)

        return self.res


def make_loops(pars, dtype='double'):
    # 0.2 ms on sparkle to form the final loops
    radius, length, thickness, axis_theta, axis_phi = \
       [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
         for base in GpuCoreShellCylinder.PD_PARS]
    parts = [
        radius.get_weights(pars['radius'], 0, 10000, True),
        length.get_weights(pars['length'], 0, 10000, True),
        thickness.get_weights(pars['thickness'], 0, 10000, True),
        axis_theta.get_weights(pars['axis_theta'], -np.inf, np.inf, False),
        axis_phi.get_weights(pars['axis_phi'], -np.inf, np.inf, False),
        ]
    # Make sure that weights are normalized to peaks at 1 so that
    # the tolerance term can be used properly on truncated distributions
    loops = np.hstack((v,w/w.max()) for v,w in parts)
    #loops = np.hstack(parts)
    loops = np.ascontiguousarray(loops.T, dtype).flatten()
    return loops, [len(p[0]) for p in parts]