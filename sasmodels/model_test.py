# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:43:56 2015

@author: David
"""

import unittest
import warnings
import numpy as np

from os.path import basename, dirname, join as joinpath
from glob import glob

try:
    from .kernelcl import load_model
except ImportError,exc:
    warnings.warn(str(exc))
    warnings.warn("using ctypes instead")
    from .kerneldll import load_model

def load_kernel(model, dtype='single'):   
    kernel = load_model(model, dtype=dtype)
    kernel.info['defaults'] = dict((p[0],p[2]) for p in kernel.info['parameters'])
    return kernel

def get_weights(model, pars, name):
    from . import weights
    
    relative = name in model.info['partype']['pd-rel']
    disperser = pars.get(name+"_pd_type", "gaussian")
    value = pars.get(name, model.info['defaults'][name])
    width = pars.get(name+"_pd", 0.0)
    npts = pars.get(name+"_pd_n", 30)
    nsigma = pars.get(name+"_pd_nsigma", 3.0)
    v,w = weights.get_weights(
            disperser, npts, width, nsigma,
            value, model.info['limits'][name], relative)
    return v,w/np.sum(w)

def eval_kernel(kernel, q, pars, cutoff=1e-5):
    input = kernel.make_input(q)
    finput = kernel(input)

    fixed_pars = [pars.get(name, finput.info['defaults'][name])
                  for name in finput.fixed_pars]
    pd_pars = [get_weights(finput, pars, p) for p in finput.pd_pars]
    return finput(fixed_pars, pd_pars, cutoff)

def annotate_exception(exc, msg):
    """
    Add an annotation to the current exception, which can then be forwarded
    to the caller using a bare "raise" statement to reraise the annotated
    exception.
    Example::
        >>> D = {}
        >>> try: 
        ...    print D['hello']
        ... except Exception,exc: 
        ...    annotate_exception(exc, "while accessing 'D'")
        ...    raise
        Traceback (most recent call last):
            ...
        KeyError: "hello while accessing 'D'"
    """
    args = exc.args
    if not args:
        arg0 = msg
    else:
        arg0 = " ".join((args[0],msg))
    exc.args = tuple([arg0] + list(args[1:]))
    
def suite():
    root = dirname(__file__)
    files = sorted(glob(joinpath(root, 'models', "[a-zA-Z]*.py")))
    models_names = [basename(f)[:-3] for f in files]
    
    suite = unittest.TestSuite()
    
    for model_name in models_names:
        module = __import__('sasmodels.models.' + model_name)
        module = getattr(module, 'models', None)

        model = getattr(module, model_name, None)
        tests = getattr(model, 'tests', [])
        
        if tests:
            #print '------'
            #print 'found tests in', model_name
            #print '------'
    
            kernel = load_kernel(model)
            suite.addTest(ModelTestCase(model_name, kernel, tests))

    return suite

class ModelTestCase(unittest.TestCase):
    
    def __init__(self, model_name, kernel, tests):
        unittest.TestCase.__init__(self)
        
        self.model_name = model_name
        self.kernel = kernel
        self.tests = tests

    def runTest(self):
        #print '------'
        #print self.model_name
        #print '------'
        try:
            for test in self.tests:
                params = test[0]
                Q = test[1]
                I = test[2]
                      
                if not isinstance(Q, list):
                    Q = [Q]
                if not isinstance(I, list):
                    I = [I]
                    
                if isinstance(Q[0], tuple):
                    npQ = [np.array([Qi[d] for Qi in Q]) for d in xrange(len(Q[0]))]
                else:
                    npQ = [np.array(Q)]

                self.assertTrue(Q)
                self.assertEqual(len(I), len(Q))    
            
                Iq = eval_kernel(self.kernel, npQ, params)
            
                self.assertGreater(len(Iq), 0)    
                self.assertEqual(len(I), len(Iq))              
                
                for q, i, iq in zip(Q, I, Iq):
                    err = np.abs(i - iq)
                    nrm = np.abs(i)
            
                    self.assertLess(err * 10**5, nrm, 'q:%s; expected:%s; actual:%s' % (q, i, iq))
                    
        except Exception,exc: 
            annotate_exception(exc, '\r\nModel: %s' % self.model_name)
            raise

def main():
    #unittest.main()
    runner = unittest.TextTestRunner()
    runner.run(suite())

if __name__ == "__main__":
    main()
