# -*- coding: utf-8 -*-
"""
Run model unit tests.

Usage::

    python -m sasmodels.model_test [opencl|dll|opencl_and_dll] model1 model2 ...

    if model1 is 'all', then all except the remaining models will be tested

Each model is tested using the default parameters at q=0.1, (qx,qy)=(0.1,0.1),
and the ER and VR are computed.  The return values at these points are not
considered.  The test is only to verify that the models run to completion,
and do not produce inf or NaN.

Tests are defined with the *tests* attribute in the model.py file.  *tests*
is a list of individual tests to run, where each test consists of the
parameter values for the test, the q-values and the expected results.  For
the effective radius test, the q-value should be 'ER'.  For the VR test,
the q-value should be 'VR'.  For 1-D tests, either specify the q value or
a list of q-values, and the corresponding I(q) value, or list of I(q) values.

That is::

    tests = [
        [ {parameters}, q, I(q)],
        [ {parameters}, [q], [I(q)] ],
        [ {parameters}, [q1, q2, ...], [I(q1), I(q2), ...]],

        [ {parameters}, (qx, qy), I(qx, Iqy)],
        [ {parameters}, [(qx1, qy1), (qx2, qy2), ...],
                        [I(qx1,qy1), I(qx2,qy2), ...]],

        [ {parameters}, 'ER', ER(pars) ],
        [ {parameters}, 'VR', VR(pars) ],
        ...
    ]

Parameters are *key:value* pairs, where key is one of the parameters of the
model and value is the value to use for the test.  Any parameters not given
in the parameter list will take on the default parameter value.

Precision defaults to 5 digits (relative).
"""

import sys
import unittest

import numpy as np

from .core import list_models, load_model_definition, load_model, HAVE_OPENCL
from .core import make_kernel, call_kernel, call_ER, call_VR
from .exception import annotate_exception


def make_suite(loaders, models):

    ModelTestCase = _hide_model_case_from_nosetests()
    suite = unittest.TestSuite()

    if models[0] == 'all':
        skip = models[1:]
        models = list_models()
    else:
        skip = []
    for model_name in models:
        if model_name in skip: continue
        model_definition = load_model_definition(model_name)

        smoke_tests = [
            [{},0.1,None],
            [{},(0.1,0.1),None],
            [{},'ER',None],
            [{},'VR',None],
            ]
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
            if not ispy and ('opencl' in loaders and HAVE_OPENCL):
                test_name = "Model: %s, Kernel: OpenCL"%model_name
                test_method = "test_%s_opencl" % model_name
                test = ModelTestCase(test_name, model_definition,
                                     tests, test_method,
                                     platform="ocl", dtype='single')
                #print "defining", test_name
                suite.addTest(test)

            # test using dll if desired
            if ispy or 'dll' in loaders:
                test_name = "Model: %s, Kernel: dll"%model_name
                test_method = "test_%s_dll" % model_name
                test = ModelTestCase(test_name, model_definition,
                                     tests, test_method,
                                     platform="dll", dtype="double")
                suite.addTest(test)

    return suite


def _hide_model_case_from_nosetests():
    class ModelTestCase(unittest.TestCase):
        def __init__(self, test_name, definition, tests, test_method,
                     platform, dtype):
            self.test_name = test_name
            self.definition = definition
            self.tests = tests
            self.platform = platform
            self.dtype = dtype

            setattr(self, test_method, self._runTest)
            unittest.TestCase.__init__(self, test_method)

        def _runTest(self):
            try:
                model = load_model(self.definition, dtype=self.dtype,
                                   platform=self.platform)
                for test in self.tests:
                    self._run_one_test(model, test)

            except Exception,exc:
                annotate_exception(exc, self.test_name)
                raise

        def _run_one_test(self, model, test):
            pars, x, y = test

            if not isinstance(y, list):
                y = [y]
            if not isinstance(x, list):
                x = [x]

            self.assertEqual(len(y), len(x))

            if x[0] == 'ER':
                actual = [call_ER(model.info, pars)]
            elif x[0] == 'VR':
                actual = [call_VR(model.info, pars)]
            elif isinstance(x[0], tuple):
                Qx,Qy = zip(*x)
                q_vectors = [np.array(Qx), np.array(Qy)]
                kernel = make_kernel(model, q_vectors)
                actual = call_kernel(kernel, pars)
            else:
                q_vectors = [np.array(x)]
                kernel = make_kernel(model, q_vectors)
                actual = call_kernel(kernel, pars)

            self.assertGreater(len(actual), 0)
            self.assertEqual(len(y), len(actual))

            for xi, yi, actual_yi in zip(x, y, actual):
                if yi is None:
                    # smoke test --- make sure it runs and produces a value
                    self.assertTrue(np.isfinite(actual_yi),
                        'invalid f(%s): %s' % (xi, actual_yi))
                else:
                    err = abs(yi - actual_yi)
                    nrm = abs(yi)
                    self.assertLess(err * 10**5, nrm,
                        'f(%s); expected:%s; actual:%s' % (xi, yi, actual_yi))

    return ModelTestCase



def main():
    """
    Run tests given is sys.argv.

    Returns 0 if success or 1 if any tests fail.
    """
    import xmlrunner

    models = sys.argv[1:]
    if models and models[0] == 'opencl':
        if not HAVE_OPENCL:
            print >>sys.stderr, "opencl is not available"
            return 1
        loaders = ['opencl']
        models = models[1:]
    elif models and models[0] == 'dll':
        # TODO: test if compiler is available?
        loaders = ['dll']
        models = models[1:]
    elif models and models[0] == 'opencl_and_dll':
        loaders = ['opencl', 'dll']
        models = models[1:]
    else:
        loaders = ['opencl', 'dll']
    if not models:
        print >>sys.stderr, "usage: python -m sasmodels.model_test [opencl|dll|opencl_and_dll] model1 model2 ..."
        print >>sys.stderr, "if model1 is 'all', then all except the remaining models will be tested"
        return 1

    #runner = unittest.TextTestRunner()
    runner = xmlrunner.XMLTestRunner(output='logs')
    result = runner.run(make_suite(loaders, models))
    return 1 if result.failures or result.errors else 0


def model_tests():
    """
    Test runner visible to nosetests.

    Run "nosetests sasmodels" on the command line to invoke it.
    """
    tests = make_suite(['opencl','dll'],['all'])
    for test_i in tests:
        yield test_i._runTest


if __name__ == "__main__":
    sys.exit(main())
