#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import datetime
import warnings

import numpy as np

def load_model(modelname):
    from os.path import abspath, join as joinpath, dirname
    from sasmodels import gen
    modelpath = abspath(joinpath(dirname(gen.__file__), 'models',
                                 modelname+'.c'))
    return gen.make(modelpath)



def opencl_model(modelname, dtype="single"):
    from sasmodels import gpu

    source, info, _ = load_model(modelname)
    # for debugging, save source to a .cl file, edit it, and reload as model
    #open(modelname+'.cl','w').write(source)
    #source = open(modelname+'.cl','r').read()
    return gpu.GpuModel(source, info, dtype)


if sys.platform == 'darwin':
    COMPILE = "gcc-mp-4.7 -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm -lgomp"
elif os.name == 'nt':
    COMPILE = "gcc -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm"
else:
    COMPILE = "cc -shared -fPIC -std=c99 -fopenmp -O2 -Wall %s -o %s -lm"
DLL_PATH = "/tmp"


def dll_path(info):
    from os.path import join as joinpath, split as splitpath, splitext
    basename = splitext(splitpath(info['filename'])[1])[0]
    return joinpath(DLL_PATH, basename+'.so')


def dll_model(modelname):
    import os
    from sasmodels import dll

    source, info, _ = load_model(modelname)
    dllpath = dll_path(info)
    if not os.path.exists(dllpath) \
            or (os.path.getmtime(dllpath) < os.path.getmtime(info['filename'])):
        # Replace with a proper temp file
        srcfile = '/tmp/%s.c'%modelname
        open(srcfile, 'w').write(source)
        os.system(COMPILE%(srcfile, dllpath))
    return dll.DllModel(dllpath, info)


TIC = None
def tic():
    global TIC
    then = datetime.datetime.now()
    TIC = lambda: (datetime.datetime.now()-then).total_seconds()
    return TIC


def toc():
    return TIC()


def load_data(filename):
    from sans.dataloader.loader import Loader
    loader = Loader()
    data = loader.load(filename)
    if data is None:
        raise IOError("Data %r could not be loaded"%filename)
    return data


def fake_data2D(qx, qy=None):
    from sans.dataloader.data_info import Data2D, Detector

    if qy is None:
        qy = qx
    Qx,Qy = np.meshgrid(qx,qy)
    Qx,Qy = Qx.flatten(), Qy.flatten()
    Iq = 100*np.ones_like(Qx)
    dIq = np.sqrt(Iq)
    mask = np.ones(len(Iq), dtype='bool')

    data = Data2D()
    data.filename = "fake data"
    data.qx_data = Qx
    data.qy_data = Qy
    data.data = Iq
    data.err_data = dIq
    data.mask = mask

    # 5% dQ/Q resolution
    data.dqx_data = 0.05*Qx
    data.dqy_data = 0.05*Qy

    detector = Detector()
    detector.pixel_size.x = 5 # mm
    detector.pixel_size.y = 5 # mm
    detector.distance = 4 # m
    data.detector.append(detector)
    data.xbins = qx
    data.ybins = qy
    data.source.wavelength = 5 # angstroms
    data.source.wavelength_unit = "A"
    data.Q_unit = "1/A"
    data.I_unit = "1/cm"
    data.q_data = np.sqrt(Qx**2 + Qy**2)
    data.xaxis("Q_x", "A^{-1}")
    data.yaxis("Q_y", "A^{-1}")
    data.zaxis("Intensity", r"\text{cm}^{-1}")
    return data


def set_beam_stop(data, radius, outer=None):
    from sans.dataloader.manipulations import Ringcut
    if hasattr(data, 'qx_data'):
        data.mask = Ringcut(0, radius)(data)
        if outer is not None:
            data.mask += Ringcut(outer,np.inf)(data)
    else:
        data.mask = (data.x>=radius)
        if outer is not None:
            data.mask &= (data.x<outer)


def set_half(data, half):
    from sans.dataloader.manipulations import Boxcut
    if half == 'right':
        data.mask += Boxcut(x_min=-np.inf, x_max=0.0, y_min=-np.inf, y_max=np.inf)(data)
    if half == 'left':
        data.mask += Boxcut(x_min=0.0, x_max=np.inf, y_min=-np.inf, y_max=np.inf)(data)


def set_top(data, max):
    from sans.dataloader.manipulations import Boxcut
    data.mask += Boxcut(x_min=-np.inf, x_max=np.inf, y_min=-np.inf, y_max=max)(data)


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
    print np.exp(np.mean(mtheory)), np.std(mtheory),np.max(mtheory),np.min(mtheory)

    plt.subplot(1, 3, 1)
    plot_data(data, mdata, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(1, 3, 2)
    plot_data(data, mtheory, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.subplot(1, 3, 3)
    print abs(mresid).max()
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

    plt.subplot(121)
    plt.errorbar(data.x, mdata, yerr=data.dy)
    plt.plot(data.x, mtheory, '-', hold=True)
    plt.yscale(view)
    plt.subplot(122)
    plt.plot(data.x, mresid, 'x')
    #plt.axhline(1, color='black', ls='--',lw=1, hold=True)
    #plt.axhline(0, color='black', lw=1, hold=True)
    #plt.axhline(-1, color='black', ls='--',lw=1, hold=True)


class BumpsModel(object):
    def __init__(self, data, model, cutoff=1e-5, **kw):
        from bumps.names import Parameter
        from . import gpu

        # interpret data
        self.is2D = hasattr(data,'qx_data')
        self.data = data
        if self.is2D:
            self.index = (data.mask==0) & (~np.isnan(data.data))
            self.iq = data.data[self.index]
            self.diq = data.err_data[self.index]
            self._theory = np.zeros_like(data.data)
            q_vectors = [data.qx_data, data.qy_data]
        else:
            self.index = (data.mask==0) & (~np.isnan(data.y))
            self.iq = data.y[self.index]
            self.diq = data.dy[self.index]
            self._theory = np.zeros_like(data.y)
            q_vectors = [data.x]
        #input = model.make_input(q_vectors)
        input = model.make_input([v[self.index] for v in q_vectors])

        # create model
        self.fn = model(input)
        self.cutoff = cutoff

        # define bumps parameters
        pars = []
        extras = []
        for p in model.info['parameters']:
            name, default, limits, ptype = p[0], p[2], p[3], p[4]
            value = kw.pop(name, default)
            setattr(self, name, Parameter.default(value, name=name, limits=limits))
            pars.append(name)
        for name in model.info['partype']['pd-2d']:
            for xpart,xdefault,xlimits in [
                    ('_pd', 0, limits),
                    ('_pd_n', 35, (0,1000)),
                    ('_pd_nsigma', 3, (0, 10)),
                    ('_pd_type', 'gaussian', None),
                ]:
                xname = name+xpart
                xvalue = kw.pop(xname, xdefault)
                if xlimits is not None:
                    xvalue = Parameter.default(xvalue, name=xname, limits=xlimits)
                    pars.append(xname)
                setattr(self, xname, xvalue)
        self._parameter_names = pars
        if kw:
            raise TypeError("unexpected parameters: %s"%(", ".join(sorted(kw.keys()))))
        self.update()

    def update(self):
        self._cache = {}

    def numpoints(self):
        return len(self.iq)

    def parameters(self):
        return dict((k,getattr(self,k)) for k in self._parameter_names)

    def theory(self):
        if 'theory' not in self._cache:
            pars = [getattr(self,p).value for p in self.fn.fixed_pars]
            pd_pars = [self._get_weights(p) for p in self.fn.pd_pars]
            #print pars
            self._theory[self.index] = self.fn(pars, pd_pars, self.cutoff)
            #self._theory[:] = self.fn.eval(pars, pd_pars)
            self._cache['theory'] = self._theory
        return self._cache['theory']

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

    def _get_weights(self, par):
        from . import weights

        relative = self.fn.info['partype']['pd-rel']
        limits = self.fn.info['limits']
        disperser,value,npts,width,nsigma = [getattr(self, par+ext)
                for ext in ('_pd_type','','_pd_n','_pd','_pd_nsigma')]
        v,w = weights.get_weights(
            disperser, int(npts.value), width.value, nsigma.value,
            value.value, limits[par], par in relative)
        return v,w/w.max()


def demo():
    data = load_data('DEC07086.DAT')
    set_beam_stop(data, 0.004)
    plot_data(data)
    import matplotlib.pyplot as plt; plt.show()


if __name__ == "__main__":
    demo()
