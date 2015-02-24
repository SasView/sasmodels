# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 11:43:56 2015

@author: David
"""

import sys
import unittest

import numpy as np

from .core import list_models, load_model_definition
from .core import load_model_cl, load_model_dll
from .core import make_kernel, call_kernel

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
        exc.args = (msg,)
    else:
        try:
            arg0 = " ".join((args[0],msg))
            exc.args = tuple([arg0] + list(args[1:]))
        except:
            exc.args = (" ".join((str(exc),msg)),)
    
def suite(loaders, models):

    suite = unittest.TestSuite()

    if models[0] == 'all':
        skip = models[1:]
        models = list_models()
    else:
        skip = []
    for model_name in models:
        if model_name in skip: continue
        model_definition = load_model_definition(model_name)

        smoke_tests = [[{},0.1,None],[{},(0.1,0.1),None]]
        tests = smoke_tests + getattr(model_definition, 'tests', [])
        
        if tests: # in case there are no smoke tests...
            #print '------'
            #print 'found tests in', model_name
            #print '------'

            # if ispy then use the dll loader to call pykernel
            # don't try to call cl kernel since it will not be
            # available in some environmentes.
            ispy = callable(getattr(model_definition,'Iq', None))

            # test using opencl if desired
            if not ispy and ('opencl' in loaders and load_model_cl):
                test_name = "Model: %s, Kernel: OpenCL"%model_name
                test = ModelTestCase(test_name, model_definition,
                                     load_model_cl, tests)
                #print "defining", test_name
                suite.addTest(test)

            # test using dll if desired
            if ispy or ('dll' in loaders and load_model_dll):
                test_name = "Model: %s, Kernel: dll"%model_name
                test = ModelTestCase(test_name, model_definition,
                                     load_model_dll, tests)
                #print "defining", test_name
                suite.addTest(test)

    return suite

class ModelTestCase(unittest.TestCase):
    
    def __init__(self, test_name, definition, loader, tests):
        unittest.TestCase.__init__(self)
        
        self.test_name = test_name
        self.definition = definition
        self.loader = loader
        self.tests = tests

    def runTest(self):
        #print "running", self.test_name
        try:
            model = self.loader(self.definition)
            for test in self.tests:
                pars, Q, I = test

                if not isinstance(Q, list):
                    Q = [Q]
                if not isinstance(I, list):
                    I = [I]
                    
                if isinstance(Q[0], tuple):
                    Qx,Qy = zip(*Q)
                    Q_vectors = [np.array(Qx), np.array(Qy)]
                else:
                    Q_vectors = [np.array(Q)]

                self.assertEqual(len(I), len(Q))

                kernel = make_kernel(model, Q_vectors)
                Iq = call_kernel(kernel, pars)
            
                self.assertGreater(len(Iq), 0)    
                self.assertEqual(len(I), len(Iq))              
                
                for q, i, iq in zip(Q, I, Iq):
                    if i is None: continue # smoke test --- make sure it runs
                    err = abs(i - iq)
                    nrm = abs(i)
            
                    self.assertLess(err * 10**5, nrm, 'q:%s; expected:%s; actual:%s' % (q, i, iq))
                    
        except Exception,exc: 
            annotate_exception(exc, self.test_name)
            raise

def main():

    models = sys.argv[1:]
    if models and models[0] == 'opencl':
        if load_model_cl is None:
            print >>sys.stderr, "opencl is not available"
            sys.exit(1)
        loaders = ['opencl']
        models = models[1:]
    elif models and models[0] == 'dll':
        # TODO: test if compiler is available?
        loaders = ['dll']
        models = models[1:]
    else:
        loaders = ['opencl', 'dll']
    if models:
        runner = unittest.TextTestRunner()
        runner.run(suite(loaders, models))
    else:
        print >>sys.stderr, "usage: python -m sasmodels.model_test [opencl|dll] model1 model2 ..."
        print >>sys.stderr, "if model1 is 'all', then all except the remaining models will be tested"

if __name__ == "__main__":
    main()
