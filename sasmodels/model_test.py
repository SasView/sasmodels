# -*- coding: utf-8 -*-
"""
Run model unit tests.

Usage::

    python -m sasmodels.model_test [opencl|dll|opencl_and_dll] model1 model2 ...

    if model1 is 'all', then all except the remaining models will be tested

Each model is tested using the default parameters at q=0.1, (qx, qy)=(0.1, 0.1),
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
                        [I(qx1, qy1), I(qx2, qy2), ...]],

        [ {parameters}, 'ER', ER(pars) ],
        [ {parameters}, 'VR', VR(pars) ],
        ...
    ]

Parameters are *key:value* pairs, where key is one of the parameters of the
model and value is the value to use for the test.  Any parameters not given
in the parameter list will take on the default parameter value.

Precision defaults to 5 digits (relative).
"""
from __future__ import print_function

import sys
import unittest

import numpy as np

from .core import list_models, load_model_info, build_model, HAVE_OPENCL
from .core import call_kernel, call_ER, call_VR
from .exception import annotate_exception

#TODO: rename to tests so that tab completion works better for models directory

def make_suite(loaders, models):
    """
    Construct the pyunit test suite.

    *loaders* is the list of kernel drivers to use, which is one of
    *["dll", "opencl"]*, *["dll"]* or *["opencl"]*.  For python models,
    the python driver is always used.

    *models* is the list of models to test, or *["all"]* to test all models.
    """

    ModelTestCase = _hide_model_case_from_nosetests()
    suite = unittest.TestSuite()

    if models[0] == 'all':
        skip = models[1:]
        models = list_models()
    else:
        skip = []
    for model_name in models:
        if model_name in skip: continue
        model_info = load_model_info(model_name)

        #print('------')
        #print('found tests in', model_name)
        #print('------')

        # if ispy then use the dll loader to call pykernel
        # don't try to call cl kernel since it will not be
        # available in some environmentes.
        is_py = callable(model_info['Iq'])

        if is_py:  # kernel implemented in python
            test_name = "Model: %s, Kernel: python"%model_name
            test_method_name = "test_%s_python" % model_name
            test = ModelTestCase(test_name, model_info,
                                 test_method_name,
                                 platform="dll",  # so that
                                 dtype="double")
            suite.addTest(test)
        else:   # kernel implemented in C
            # test using opencl if desired and available
            if 'opencl' in loaders and HAVE_OPENCL:
                test_name = "Model: %s, Kernel: OpenCL"%model_name
                test_method_name = "test_%s_opencl" % model_name
                # Using dtype=None so that the models that are only
                # correct for double precision are not tested using
                # single precision.  The choice is determined by the
                # presence of *single=False* in the model file.
                test = ModelTestCase(test_name, model_info,
                                     test_method_name,
                                     platform="ocl", dtype=None)
                #print("defining", test_name)
                suite.addTest(test)

            # test using dll if desired
            if 'dll' in loaders:
                test_name = "Model: %s, Kernel: dll"%model_name
                test_method_name = "test_%s_dll" % model_name
                test = ModelTestCase(test_name, model_info,
                                     test_method_name,
                                     platform="dll",
                                     dtype="double")
                suite.addTest(test)

    return suite


def _hide_model_case_from_nosetests():
    class ModelTestCase(unittest.TestCase):
        """
        Test suit for a particular model with a particular kernel driver.

        The test suite runs a simple smoke test to make sure the model
        functions, then runs the list of tests at the bottom of the model
        description file.
        """
        def __init__(self, test_name, model_info, test_method_name,
                     platform, dtype):
            self.test_name = test_name
            self.info = model_info
            self.platform = platform
            self.dtype = dtype

            setattr(self, test_method_name, self._runTest)
            unittest.TestCase.__init__(self, test_method_name)

        def _runTest(self):
            smoke_tests = [
                # test validity at reasonable values
                [{}, 0.1, None],
                [{}, (0.1, 0.1), None],
                # test validity at q = 0
                #[{}, 0.0, None],
                #[{}, (0.0, 0.0), None],
                # test that ER/VR will run if they exist
                [{}, 'ER', None],
                [{}, 'VR', None],
                ]

            tests = self.info['tests']
            try:
                model = build_model(self.info, dtype=self.dtype,
                                    platform=self.platform)
                for test in smoke_tests + tests:
                    self._run_one_test(model, test)

                if not tests and self.platform == "dll":
                    ## Uncomment the following to make forgetting the test
                    ## values an error.  Only do so for the "dll" tests
                    ## to reduce noise from both opencl and dll, and because
                    ## python kernels use platform="dll".
                    #raise Exception("No test cases provided")
                    pass

            except:
                annotate_exception(self.test_name)
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
                Qx, Qy = zip(*x)
                q_vectors = [np.array(Qx), np.array(Qy)]
                kernel = model.make_kernel(q_vectors)
                actual = call_kernel(kernel, pars)
            else:
                q_vectors = [np.array(x)]
                kernel = model.make_kernel(q_vectors)
                actual = call_kernel(kernel, pars)

            self.assertTrue(len(actual) > 0)
            self.assertEqual(len(y), len(actual))

            for xi, yi, actual_yi in zip(x, y, actual):
                if yi is None:
                    # smoke test --- make sure it runs and produces a value
                    self.assertTrue(np.isfinite(actual_yi),
                                    'invalid f(%s): %s' % (xi, actual_yi))
                elif np.isnan(yi):
                    self.assertTrue(np.isnan(actual_yi),
                                    'f(%s): expected:%s; actual:%s'
                                    % (xi, yi, actual_yi))
                else:
                    # is_near does not work for infinite values, so also test
                    # for exact values.  Note that this will not
                    self.assertTrue(yi==actual_yi or is_near(yi, actual_yi, 5),
                                    'f(%s); expected:%s; actual:%s'
                                    % (xi, yi, actual_yi))

    return ModelTestCase

def is_near(target, actual, digits=5):
    """
    Returns true if *actual* is within *digits* significant digits of *target*.
    """
    import math
    shift = 10**math.ceil(math.log10(abs(target)))
    return abs(target-actual)/shift < 1.5*10**-digits

def main():
    """
    Run tests given is sys.argv.

    Returns 0 if success or 1 if any tests fail.
    """
    try:
        from xmlrunner import XMLTestRunner as TestRunner
        test_args = { 'output': 'logs' }
    except ImportError:
        from unittest import TextTestRunner as TestRunner
        test_args = { }

    models = sys.argv[1:]
    if models and models[0] == '-v':
        verbosity = 2
        models = models[1:]
    else:
        verbosity = 1
    if models and models[0] == 'opencl':
        if not HAVE_OPENCL:
            print("opencl is not available")
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
        print("""\
usage:
  python -m sasmodels.model_test [-v] [opencl|dll] model1 model2 ...

If -v is included on the command line, then use verboe output.

If neither opencl nor dll is specified, then models will be tested with
both opencl and dll; the compute target is ignored for pure python models.

If model1 is 'all', then all except the remaining models will be tested.

""")

        return 1

    runner = TestRunner(verbosity=verbosity, **test_args)
    result = runner.run(make_suite(loaders, models))
    return 1 if result.failures or result.errors else 0


def model_tests():
    """
    Test runner visible to nosetests.

    Run "nosetests sasmodels" on the command line to invoke it.
    """
    tests = make_suite(['opencl', 'dll'], ['all'])
    for test_i in tests:
        yield test_i._runTest


if __name__ == "__main__":
    sys.exit(main())
