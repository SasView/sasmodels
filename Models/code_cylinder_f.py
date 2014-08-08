#!/usr/bin/env python
# -*- coding: utf-8 -*-


import ctypes
from ctypes import c_int, c_double, c_void_p
import numpy as np
import pyopencl as cl
from pyopencl import mem_flags as mf

from weights import GaussianDispersion
from sasmodel import card, set_precision, set_precision_1d, tic, toc


class GpuCylinder(object):
    PARS = {
        'scale':1,'radius':1,'length':1,'sldCyl':1e-6,'sldSolv':0,'background':0,
        'cyl_theta':0,'cyl_phi':0,
    }
    PD_PARS = ['radius', 'length', 'cyl_theta', 'cyl_phi']

    def __init__(self, qx, qy, dtype='float32', cutoff=1e-5):

        #create context, queue, and build program
        ctx,_queue = card()
        bessel = open('Kernel/NR_BessJ1.cpp').read()
        kernel = open('Kernel/Kernel-Cylinder_f.cpp').read()
        src, qx, qy = set_precision("\n".join((bessel,kernel)), qx, qy, dtype=dtype)
        self.prg = cl.Program(ctx, src).build()
        self.qx, self.qy = qx, qy
        self.cutoff = cutoff

        #buffers
        self.qx_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(ctx, cl.mem_flags.READ_WRITE, self.qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):
        tic()

        ctx,queue = card()
        real = np.float32 if self.qx.dtype == np.dtype('float32') else np.float64
        loops, loop_lengths = make_loops(pars, dtype=self.qx.dtype)
        loops_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=loops)
        loops_l = cl.LocalMemory(len(loops.data))

        self.prg.CylinderKernel(queue, self.qx.shape, None,
            self.qx_b, self.qy_b, self.res_b,
            loops_b, loops_l, real(self.cutoff),
            real(pars['scale']), real(pars['background']),
            real(pars['sldCyl']-pars['sldSolv']),
            *[np.uint32(pn) for pn in loop_lengths])

        cl.enqueue_copy(queue, self.res, self.res_b)
        print toc()*1000, self.qx.shape[0]

        return self.res

class CpuCylinder(GpuCylinder): #inherit parameters only
    def __init__(self, qx, qy, dtype='float32', cutoff=1e-5):
        self.qx, self.qy = [np.ascontiguousarray(v,'d') for v in qx,qy]
        self.cutoff = cutoff
        self.res = np.empty_like(self.qx)
        self.dll = ctypes.CDLL('Kernel/cylinder.so')
        self.fn = self.dll['CylinderKernel']
        self.fn.argtypes = [
            c_void_p, c_void_p, c_void_p, c_int,
            c_void_p, c_double, c_double, c_double, c_double,
            c_int, c_int, c_int, c_int
        ]
    def eval(self, pars):
        loops, loop_lengths = make_loops(pars, dtype=self.qx.dtype)
        self.fn(self.qx.ctypes.data, self.qy.ctypes.data, self.res.ctypes.data, len(self.qx),
            loops.ctypes.data, self.cutoff, pars['scale'], pars['background'],
            pars['sldCyl']-pars['sldSolv'], *loop_lengths)

        return self.res

def make_loops(pars, dtype='double'):
    # 0.2 ms on sparkle to form the final loops
    radius, length, theta, phi = \
       [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
         for base in GpuCylinder.PD_PARS]
    parts = [
        radius.get_weights(pars['radius'], 0, 10000, True),
        length.get_weights(pars['length'], 0, 10000, True),
        theta.get_weights(pars['cyl_theta'], -np.inf, np.inf, False),
        phi.get_weights(pars['cyl_phi'], -np.inf, np.inf, False),
        ]
    # Make sure that weights are normalized to peaks at 1 so that
    # the tolerance term can be used properly on truncated distributions
    loops = np.hstack((v,w/w.max()) for v,w in parts)
    #loops = np.hstack(parts)
    loops = np.ascontiguousarray(loops.T, dtype).flatten()
    return loops, [len(p[0]) for p in parts]

class OneDGpuCylinder(object):
    PARS = {
        'scale':1,'radius':1,'length':1,'sldCyl':1e-6,'sldSolv':0,'background':0,
        'bolim':0, 'uplim':90
    }
    PD_PARS = ['radius', 'length']

    def __init__(self, q, dtype='float32'):

        #create context, queue, and build program
        ctx,_queue = card()
        trala = open('Kernel/NR_BessJ1.cpp').read()+"\n"+open('Kernel/OneDCyl_Kfun.cpp').read()+"\n"+open('Kernel/Kernel-OneDCylinder.cpp').read()
        src, self.q = set_precision_1d(trala, q, dtype=dtype)
        self.prg = cl.Program(ctx, src).build()

        #buffers
        mf = cl.mem_flags
        self.q_b = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.q)
        self.res_b = cl.Buffer(ctx, mf.WRITE_ONLY, q.nbytes)
        self.res = np.empty_like(self.q)

    def eval(self, pars):

        _ctx,queue = card()
        radius, length = \
            [GaussianDispersion(int(pars[base+'_pd_n']), pars[base+'_pd'], pars[base+'_pd_nsigma'])
             for base in OneDGpuCylinder.PD_PARS]

        #Get the weights for each
        radius.value, radius.weight = radius.get_weights(pars['radius'], 0, 10000, True)
        length.value, length.weight = length.get_weights(pars['length'], 0, 10000, True)

        #Perform the computation, with all weight points
        sum, norm, vol = 0.0, 0.0, 0.0,
        sub = pars['sldCyl'] - pars['sldSolv']

        real = np.float32 if self.q.dtype == np.dtype('float32') else np.float64
        #Loop over radius, length, theta, phi weight points
        for r in xrange(len(radius.weight)):
            for l in xrange(len(length.weight)):
                        self.prg.OneDCylKernel(queue, self.q.shape, None, self.q_b, self.res_b, real(sub),
                                           real(length.value[l]), real(radius.value[r]), real(pars['scale']),
                                           np.uint32(self.q.size), real(pars['uplim']), real(pars['bolim']))
                        cl.enqueue_copy(queue, self.res, self.res_b)
                        sum += radius.weight[r]*length.weight[l]*self.res*pow(radius.value[r],2)*length.value[l]
                        vol += radius.weight[r]*length.weight[l] *pow(radius.value[r],2)*length.value[l]
                        norm += radius.weight[r]*length.weight[l]

        if vol != 0.0 and norm != 0.0:
            sum *= norm/vol

        return sum/norm + pars['background']
