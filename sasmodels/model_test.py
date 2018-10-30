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
import traceback

try:
    from StringIO import StringIO
except ImportError:
    # StringIO.StringIO renamed to io.StringIO in Python 3
    # Note: io.StringIO exists in python 2, but using unicode instead of str
    from io import StringIO

import numpy as np  # type: ignore

from . import core
from .core import list_models, load_model_info, build_model
from .direct_model import call_kernel, call_ER, call_VR
from .exception import annotate_exception
from .modelinfo import expand_pars
from .kernelcl import use_opencl

# pylint: disable=unused-import
try:
    from typing import List, Iterator, Callable
except ImportError:
    pass
else:
    from .modelinfo import ParameterTable, ParameterSet, TestCondition, ModelInfo
    from .kernel import KernelModel
# pylint: enable=unused-import

def make_suite(loaders, models):
    # type: (List[str], List[str]) -> unittest.TestSuite
    """
    Construct the pyunit test suite.

    *loaders* is the list of kernel drivers to use, which is one of
    *["dll", "opencl"]*, *["dll"]* or *["opencl"]*.  For python models,
    the python driver is always used.

    *models* is the list of models to test, or *["all"]* to test all models.
    """
    suite = unittest.TestSuite()

    if models[0] in core.KINDS:
        skip = models[1:]
        models = list_models(models[0])
    else:
        skip = []
    for model_name in models:
        if model_name not in skip:
            model_info = load_model_info(model_name)
            _add_model_to_suite(loaders, suite, model_info)

    return suite

def _add_model_to_suite(loaders, suite, model_info):
    ModelTestCase = _hide_model_case_from_nose()

    #print('------')
    #print('found tests in', model_name)
    #print('------')

    # if ispy then use the dll loader to call pykernel
    # don't try to call cl kernel since it will not be
    # available in some environmentes.
    is_py = callable(model_info.Iq)

    # Some OpenCL drivers seem to be flaky, and are not producing the
    # expected result.  Since we don't have known test values yet for
    # all of our models, we are instead going to compare the results
    # for the 'smoke test' (that is, evaluation at q=0.1 for the default
    # parameters just to see that the model runs to completion) between
    # the OpenCL and the DLL.  To do this, we define a 'stash' which is
    # shared between OpenCL and DLL tests.  This is just a list.  If the
    # list is empty (which it will be when DLL runs, if the DLL runs
    # first), then the results are appended to the list.  If the list
    # is not empty (which it will be when OpenCL runs second), the results
    # are compared to the results stored in the first element of the list.
    # This is a horrible stateful hack which only makes sense because the
    # test suite is thrown away after being run once.
    stash = []

    if is_py:  # kernel implemented in python
        test_name = "%s-python"%model_info.name
        test_method_name = "test_%s_python" % model_info.id
        test = ModelTestCase(test_name, model_info,
                                test_method_name,
                                platform="dll",  # so that
                                dtype="double",
                                stash=stash)
        suite.addTest(test)
    else:   # kernel implemented in C

        # test using dll if desired
        if 'dll' in loaders or not use_opencl():
            test_name = "%s-dll"%model_info.name
            test_method_name = "test_%s_dll" % model_info.id
            test = ModelTestCase(test_name, model_info,
                                    test_method_name,
                                    platform="dll",
                                    dtype="double",
                                    stash=stash)
            suite.addTest(test)

        # test using opencl if desired and available
        if 'opencl' in loaders and use_opencl():
            test_name = "%s-opencl"%model_info.name
            test_method_name = "test_%s_opencl" % model_info.id
            # Using dtype=None so that the models that are only
            # correct for double precision are not tested using
            # single precision.  The choice is determined by the
            # presence of *single=False* in the model file.
            test = ModelTestCase(test_name, model_info,
                                    test_method_name,
                                    platform="ocl", dtype=None,
                                    stash=stash)
            #print("defining", test_name)
            suite.addTest(test)


def _hide_model_case_from_nose():
    # type: () -> type
    class ModelTestCase(unittest.TestCase):
        """
        Test suit for a particular model with a particular kernel driver.

        The test suite runs a simple smoke test to make sure the model
        functions, then runs the list of tests at the bottom of the model
        description file.
        """
        def __init__(self, test_name, model_info, test_method_name,
                     platform, dtype, stash):
            # type: (str, ModelInfo, str, str, DType, List[Any]) -> None
            self.test_name = test_name
            self.info = model_info
            self.platform = platform
            self.dtype = dtype
            self.stash = stash  # container for the results of the first run

            setattr(self, test_method_name, self.run_all)
            unittest.TestCase.__init__(self, test_method_name)

        def run_all(self):
            # type: () -> None
            """
            Run all the tests in the test suite, including smoke tests.
            """
            smoke_tests = [
                # test validity at reasonable values
                ({}, 0.1, None),
                ({}, (0.1, 0.1), None),
                # test validity at q = 0
                #({}, 0.0, None),
                #({}, (0.0, 0.0), None),
                # test vector form
                ({}, [0.001, 0.01, 0.1], [None]*3),
                ({}, [(0.1, 0.1)]*2, [None]*2),
                # test that ER/VR will run if they exist
                ({}, 'ER', None),
                ({}, 'VR', None),
                ]
            tests = smoke_tests
            #tests = []
            if self.info.tests is not None:
                tests += self.info.tests
            try:
                model = build_model(self.info, dtype=self.dtype,
                                    platform=self.platform)
                results = [self.run_one(model, test) for test in tests]
                if self.stash:
                    for test, target, actual in zip(tests, self.stash[0], results):
                        assert np.all(abs(target-actual) < 5e-5*abs(actual)), \
                            ("GPU/CPU comparison expected %s but got %s for %s"
                             % (target, actual, test[0]))
                else:
                    self.stash.append(results)

                # Check for missing tests.  Only do so for the "dll" tests
                # to reduce noise from both opencl and dll, and because
                # python kernels use platform="dll".
                if self.platform == "dll":
                    missing = []
                    ## Uncomment the following to require test cases
                    #missing = self._find_missing_tests()
                    if missing:
                        raise ValueError("Missing tests for "+", ".join(missing))

            except:
                annotate_exception(self.test_name)
                raise

        def _find_missing_tests(self):
            # type: () -> None
            """make sure there are 1D, 2D, ER and VR tests as appropriate"""
            model_has_VR = callable(self.info.VR)
            model_has_ER = callable(self.info.ER)
            model_has_1D = True
            model_has_2D = any(p.type == 'orientation'
                               for p in self.info.parameters.kernel_parameters)

            # Lists of tests that have a result that is not None
            single = [test for test in self.info.tests
                      if not isinstance(test[2], list) and test[2] is not None]
            tests_has_VR = any(test[1] == 'VR' for test in single)
            tests_has_ER = any(test[1] == 'ER' for test in single)
            tests_has_1D_single = any(isinstance(test[1], float) for test in single)
            tests_has_2D_single = any(isinstance(test[1], tuple) for test in single)

            multiple = [test for test in self.info.tests
                        if isinstance(test[2], list)
                        and not all(result is None for result in test[2])]
            tests_has_1D_multiple = any(isinstance(test[1][0], float)
                                        for test in multiple)
            tests_has_2D_multiple = any(isinstance(test[1][0], tuple)
                                        for test in multiple)

            missing = []
            if model_has_VR and not tests_has_VR:
                missing.append("VR")
            if model_has_ER and not tests_has_ER:
                missing.append("ER")
            if model_has_1D and not (tests_has_1D_single or tests_has_1D_multiple):
                missing.append("1D")
            if model_has_2D and not (tests_has_2D_single or tests_has_2D_multiple):
                missing.append("2D")

            return missing

        def run_one(self, model, test):
            # type: (KernelModel, TestCondition) -> None
            """Run a single test case."""
            user_pars, x, y = test
            pars = expand_pars(self.info.parameters, user_pars)
            invalid = invalid_pars(self.info.parameters, pars)
            if invalid:
                raise ValueError("Unknown parameters in test: " + ", ".join(invalid))

            if not isinstance(y, list):
                y = [y]
            if not isinstance(x, list):
                x = [x]

            self.assertEqual(len(y), len(x))

            if x[0] == 'ER':
                actual = np.array([call_ER(model.info, pars)])
            elif x[0] == 'VR':
                actual = np.array([call_VR(model.info, pars)])
            elif isinstance(x[0], tuple):
                qx, qy = zip(*x)
                q_vectors = [np.array(qx), np.array(qy)]
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
                    self.assertTrue(not np.isnan(actual_yi),
                                    'invalid f(%s): %s' % (xi, actual_yi))
                elif np.isnan(yi):
                    self.assertTrue(np.isnan(actual_yi),
                                    'f(%s): expected:%s; actual:%s'
                                    % (xi, yi, actual_yi))
                else:
                    # is_near does not work for infinite values, so also test
                    # for exact values.  Note that this will not
                    self.assertTrue(yi == actual_yi or is_near(yi, actual_yi, 5),
                                    'f(%s); expected:%s; actual:%s'
                                    % (xi, yi, actual_yi))
            return actual

    return ModelTestCase

def invalid_pars(partable, pars):
    # type: (ParameterTable, Dict[str, float])
    """
    Return a list of parameter names that are not part of the model.
    """
    names = set(p.id for p in partable.call_parameters)
    invalid = []
    for par in sorted(pars.keys()):
        parts = par.split('_pd')
        if len(parts) > 1 and parts[1] not in ("", "_n", "nsigma", "type"):
            invalid.append(par)
            continue
        if parts[0] not in names:
            invalid.append(par)
    return invalid


def is_near(target, actual, digits=5):
    # type: (float, float, int) -> bool
    """
    Returns true if *actual* is within *digits* significant digits of *target*.
    """
    import math
    shift = 10**math.ceil(math.log10(abs(target)))
    return abs(target-actual)/shift < 1.5*10**-digits

# CRUFT: old interface; should be deprecated and removed
def run_one(model_name):
    # msg = "use check_model(model_info) rather than run_one(model_name)"
    # warnings.warn(msg, category=DeprecationWarning, stacklevel=2)
    try:
        model_info = load_model_info(model_name)
    except Exception:
        output = traceback.format_exc()
        return output

    success, output = check_model(model_info)
    return output

def check_model(model_info):
    # type: (ModelInfo) -> str
    """
    Run the tests for a single model, capturing the output.

    Returns success status and the output string.
    """
    # Note that running main() directly did not work from within the
    # wxPython pycrust console.  Instead of the results appearing in the
    # window they were printed to the underlying console.
    from unittest.runner import TextTestResult, _WritelnDecorator

    # Build a object to capture and print the test results
    stream = _WritelnDecorator(StringIO())  # Add writeln() method to stream
    verbosity = 2
    descriptions = True
    result = TextTestResult(stream, descriptions, verbosity)

    # Build a test suite containing just the model
    loaders = ['opencl'] if use_opencl() else ['dll']
    suite = unittest.TestSuite()
    _add_model_to_suite(loaders, suite, model_info)

    # Warn if there are no user defined tests.
    # Note: the test suite constructed above only has one test in it, which
    # runs through some smoke tests to make sure the model runs, then runs
    # through the input-output pairs given in the model definition file.  To
    # check if any such pairs are defined, therefore, we just need to check if
    # they are in the first test of the test suite.  We do this with an
    # iterator since we don't have direct access to the list of tests in the
    # test suite.
    # In Qt5 suite.run() will clear all tests in the suite after running
    # with no way of retaining them for the test below, so let's check
    # for user tests before running the suite.
    for test in suite:
        if not test.info.tests:
            stream.writeln("Note: %s has no user defined tests."%model_info.name)
        break
    else:
        stream.writeln("Note: no test suite created --- this should never happen")

    # Run the test suite
    suite.run(result)

    # Print the failures and errors
    for _, tb in result.errors:
        stream.writeln(tb)
    for _, tb in result.failures:
        stream.writeln(tb)

    output = stream.getvalue()
    stream.close()
    return result.wasSuccessful(), output


def main(*models):
    # type: (*str) -> int
    """
    Run tests given is models.

    Returns 0 if success or 1 if any tests fail.
    """
    try:
        from xmlrunner import XMLTestRunner as TestRunner
        test_args = {'output': 'logs'}
    except ImportError:
        from unittest import TextTestRunner as TestRunner
        test_args = {}

    if models and models[0] == '-v':
        verbosity = 2
        models = models[1:]
    else:
        verbosity = 1
    if models and models[0] == 'opencl':
        if not use_opencl():
            print("opencl is not available")
            return 1
        loaders = ['opencl']
        models = models[1:]
    elif models and models[0] == 'dll':
        # TODO: test if compiler is available?
        loaders = ['dll']
        models = models[1:]
    elif models and models[0] == 'opencl_and_dll':
        loaders = ['opencl', 'dll'] if use_opencl() else ['dll']
        models = models[1:]
    else:
        loaders = ['opencl', 'dll'] if use_opencl() else ['dll']
    if not models:
        print("""\
usage:
  python -m sasmodels.model_test [-v] [opencl|dll] model1 model2 ...

If -v is included on the command line, then use verbose output.

If neither opencl nor dll is specified, then models will be tested with
both OpenCL and dll; the compute target is ignored for pure python models.

If model1 is 'all', then all except the remaining models will be tested.

""")

        return 1

    runner = TestRunner(verbosity=verbosity, **test_args)
    result = runner.run(make_suite(loaders, models))
    return 1 if result.failures or result.errors else 0


def model_tests():
    # type: () -> Iterator[Callable[[], None]]
    """
    Test runner visible to nosetests.

    Run "nosetests sasmodels" on the command line to invoke it.
    """
    loaders = ['opencl', 'dll'] if use_opencl() else ['dll']
    tests = make_suite(loaders, ['all'])
    def build_test(test):
        # In order for nosetest to show the test name, wrap the test.run_all
        # instance in function that takes the test name as a parameter which
        # will be displayed when the test is run.  Do this as a function so
        # that it properly captures the context for tests that captured and
        # run later.  If done directly in the for loop, then the looping
        # variable test will be shared amongst all the tests, and we will be
        # repeatedly testing vesicle.

        # Note: in sasview sas.sasgui.perspectives.fitting.gpu_options
        # requires that the test.description field be set.
        wrap = lambda: test.run_all()
        wrap.description = test.test_name
        return wrap
        # The following would work with nosetests and pytest:
        #     return lambda name: test.run_all(), test.test_name

    for test in tests:
        yield build_test(test)


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
