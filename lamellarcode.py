#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import pyopencl as cl
from weights import GaussianDispersion


class GpuLamellar(object):
    PARS = {
        'scale':1, 'bi_thick':1, 'sld_bi':1e-6, 'sld_sol':0, 'background':0,
    }
    PD_PARS = ['bi_thick']

    def __init__(self, qx, qy):

        self.qx = np.asarray(qx, np.float32)
        self.qy = np.asarray(qy, np.float32)
        #create context, queue, and build program
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        self.prg = cl.Program(self.ctx, open('Kernel-Lamellar.cpp').read()).build()

        #buffers
        mf = cl.mem_flags
        self.qx_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qx)
        self.qy_b = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.qy)
        self.res_b = cl.Buffer(self.ctx, mf.WRITE_ONLY, qx.nbytes)
        self.res = np.empty_like(self.qx)

    def eval(self, pars):

        bi_thick = GaussianDispersion(int(pars['bi_thick_pd_n']), pars['bi_thick_pd'], pars['bi_thick_pd_nsigma'])
        bi_thick.value, bi_thick.weight = bi_thick.get_weights(pars['bi_thick'], 0, 1000, True)

        sum, norm = 0.0, 0.0
        sub = pars['sld_bi'] - pars['sld_sol']

        for i in xrange(len(bi_thick.weight)):
            self.prg.LamellarKernel(self.queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b, np.float32(bi_thick.value[i]),
                                    np.float32(pars['scale']), np.float32(sub), np.float32(pars['background']), np.uint32(self.qx.size))
            cl.enqueue_copy(self.queue, self.res, self.res_b)

            sum += bi_thick.weight[i]*self.res
            norm += bi_thick.weight[i]

        return sum/norm + pars['background']

    def lamellar_fit(self, pars, b_n=10, b_w=.1, sigma=3):

        bi_thick = GaussianDispersion(b_n, b_w, sigma)
        bi_thick.value, bi_thick.weight = bi_thick.get_weights(pars.bi_thick, 0, 1000, True)

        sum, norm = 0.0, 0.0

        for i in xrange(len(bi_thick.weight)):
            self.prg.LamellarKernel(self.queue, self.qx.shape, None, self.qx_b, self.qy_b, self.res_b, np.float32(bi_thick.value[i]),
                                    np.float32(pars.scale), np.float32(pars.sld_bi), np.float32(pars.sld_sol),
                                    np.float32(pars.background), np.uint32(self.qx.size))
            cl.enqueue_copy(self.queue, self.res, self.res_b)

            sum += bi_thick.weight[i]*self.res
            norm += bi_thick.weight[i]

        return sum/norm + pars.background


def demo():

    from time import time
    import matplotlib.pyplot as plt

    #create qx and qy evenly spaces
    qx = np.linspace(-.01, .01, 128)
    qy = np.linspace(-.01, .01, 128)
    qx, qy = np.meshgrid(qx, qy)

    #saved shape of qx
    r_shape = qx.shape

    #reshape for calculation; resize as float32
    qx = qx.flatten()
    qy = qy.flatten()

    pars = LamellarParameters(scale=1, bi_thick=100, sld_bi=.291e-6, sld_sol=5.77e-6, background=0)
    t = time()
    result = GpuLamellar(qx, qy)
    result.x = result.lamellar_fit(pars, b_n=35, b_w=.1, sigma=3)
    result.x = np.reshape(result.x, r_shape)
    tt = time()

    print("Time taken: %f" % (tt - t))

    f = open("r.txt", "w")
    for x in xrange(len(r_shape)):
        f.write(str(result.x[x]))
        f.write("\n")

    plt.pcolormesh(np.log10(result.x))
    plt.show()


if __name__ == "__main__":
    demo()
































