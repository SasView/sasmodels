#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pyopencl as cl
from bumps.names import Parameter
from sans.dataloader.loader import Loader
from sans.dataloader.manipulations import Ringcut, Boxcut


def load_data(filename):
    loader = Loader()
    data = loader.load(filename)
    if data is None:
        raise IOError("Data %r could not be loaded"%filename)
    return data


def set_beam_stop(data, radius, outer=None):
    if hasattr(data, 'qx_data'):
        data.mask = Ringcut(0, radius)(data)
        if outer is not None:
            data.mask += Ringcut(outer,np.inf)(data)
    else:
        data.mask = (data.x>=radius)
        if outer is not None:
            data.mask &= (data.x<outer)

def set_half(data, half):
    if half == 'left':
        data.mask += Boxcut(x_min=-np.inf, x_max=0.0, y_min=-np.inf, y_max=np.inf)(data)
    if half == 'right':
        data.mask += Boxcut(x_min=0.0, x_max=np.inf, y_min=-np.inf, y_max=np.inf)(data)



def plot_data(data, iq, vmin=None, vmax=None):
    from numpy.ma import masked_array
    import matplotlib.pyplot as plt
    img = masked_array(iq, data.mask)
    xmin, xmax = min(data.qx_data), max(data.qx_data)
    ymin, ymax = min(data.qy_data), max(data.qy_data)
    plt.imshow(img.reshape(128,128),
               interpolation='nearest', aspect=1, origin='upper',
               extent=[xmin, xmax, ymin, ymax], vmin=vmin, vmax=vmax)

def plot_result2D(data, theory, view='linear'):
    import matplotlib.pyplot as plt
    from numpy.ma import masked_array, masked
    #print "not a number",sum(np.isnan(data.data))
    #data.data[data.data<0.05] = 0.5
    mdata = masked_array(data.data, data.mask)
    mdata[np.isnan(mdata)] = masked
    if view is 'log':
        mdata[mdata <= 0] = masked
        mdata = np.log10(mdata)
        mtheory = masked_array(np.log10(theory), mdata.mask)
    else:
        mtheory = masked_array(theory, mdata.mask)
    mresid = masked_array((theory-data.data)/data.err_data, data.mask)
    vmin = min(mdata.min(), mtheory.min())
    vmax = max(mdata.max(), mtheory.max())

    plt.subplot(1, 3, 1)
    plot_data(data, mdata, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(1, 3, 2)
    plot_data(data, mtheory, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(1, 3, 3)
    plot_data(data, mresid)
    plt.colorbar()


def plot_result1D(data, theory, view='linear'):
    import matplotlib.pyplot as plt
    from numpy.ma import masked_array, masked
    #print "not a number",sum(np.isnan(data.y))
    #data.y[data.y<0.05] = 0.5
    mdata = masked_array(data.y, data.mask)
    mdata[np.isnan(mdata)] = masked
    if view is 'log':
        mdata[mdata <= 0] = masked
    mtheory = masked_array(theory, mdata.mask)
    mresid = masked_array((theory-data.y)/data.dy, mdata.mask)

    plt.subplot(1,2,1)
    plt.errorbar(data.x, mdata, yerr=data.dy)
    plt.plot(data.x, mtheory, '-', hold=True)
    plt.yscale(view)
    plt.subplot(1, 2, 2)
    plt.plot(data.x, mresid, 'x')
    #plt.axhline(1, color='black', ls='--',lw=1, hold=True)
    #plt.axhline(0, color='black', lw=1, hold=True)
    #plt.axhline(-1, color='black', ls='--',lw=1, hold=True)



GPU_CONTEXT = None
GPU_QUEUE = None
def card():
    global GPU_CONTEXT, GPU_QUEUE
    if GPU_CONTEXT is None:
        GPU_CONTEXT = cl.create_some_context()
        GPU_QUEUE = cl.CommandQueue(GPU_CONTEXT)
    return GPU_CONTEXT, GPU_QUEUE


class SasModel(object):
    def __init__(self, data, model, dtype='float32', **kw):
        self.__dict__['_parameters'] = {}
        self.is2D = hasattr(data,'qx_data')
        self.data = data
        if self.is2D:
            self.index = (data.mask==0) & (~np.isnan(data.data))
            self.iq = data.data[self.index]
            self.diq = data.err_data[self.index]
            self.qx = data.qx_data
            self.qy = data.qy_data
            self.gpu = model(self.qx, self.qy, dtype=dtype)
        else:
            self.index = (data.mask==0) & (~np.isnan(data.y))
            self.iq = data.y[self.index]
            self.diq = data.dy[self.index]
            self.q = data.x
            self.gpu = model(self.q, dtype=dtype)
        pd_pars = set(base+attr for base in model.PD_PARS for attr in ('_pd','_pd_n','_pd_nsigma'))
        total_pars = set(model.PARS.keys()) | pd_pars
        extra_pars = set(kw.keys()) - total_pars
        if extra_pars:
            raise TypeError("unexpected parameters %s"%(str(extra_pars,)))
        pars = model.PARS.copy()
        pars.update((base+'_pd', 0) for base in model.PD_PARS)
        pars.update((base+'_pd_n', 35) for base in model.PD_PARS)
        pars.update((base+'_pd_nsigma', 3) for base in model.PD_PARS)
        pars.update(kw)
        self._parameters = dict((k, Parameter.default(v, name=k)) for k, v in pars.items())

    def numpoints(self):
        return len(self.iq)

    def parameters(self):
        return self._parameters

    def __getattr__(self, par):
        return self._parameters[par]

    def __setattr__(self, par, val):
        if par in self._parameters:
            self._parameters[par] = val
        else:
            self.__dict__[par] = val

    def theory(self):
        pars = dict((k,v.value) for k,v in self._parameters.items())
        result = self.gpu.eval(pars)
        return result

    def residuals(self):
        #if np.any(self.err ==0): print "zeros in err"
        return (self.theory()[self.index]-self.iq)/self.diq

    def nllf(self):
        R = self.residuals()
        #if np.any(np.isnan(R)): print "NaN in residuals"
        return 0.5*np.sum(R**2)

    def __call__(self):
        return 2*self.nllf()/self.dof

    def plot(self, view='log'):
        if self.is2D:
            plot_result2D(self.data, self.theory(), view=view)
        else:
            plot_result1D(self.data, self.theory(), view=view)

    def save(self, basename):
        pass

    def update(self):
        pass


def demo():
    data = load_data('DEC07086.DAT')
    set_beam_stop(data, 0.004)
    plot_data(data)
    import matplotlib.pyplot as plt; plt.show()

if __name__ == "__main__":
    demo()
