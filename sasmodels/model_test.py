# -*- coding: utf-8 -*-
"""
Run model unit tests.

Usage::

    python -m sasmodels.model_test [opencl|cuda|dll|all] model1 model2 ...

If model1 is 'all', then all except the remaining models will be tested.
Subgroups are also possible, such as 'py', 'single' or '1d'.  See
:func:`.core.list_models` for details.

Each model is tested using the default parameters at q=0.1, (qx, qy)=(0.1, 0.1),
and Fq is called to make sure R_eff, volume and volume ratio are computed.
The return values at these points are not considered.  The test is only to
verify that the models run to completion, and do not produce inf or NaN.

Tests are defined with the *tests* attribute in the model.py file.  *tests*
is a list of individual tests to run, where each test consists of the
parameter values for the test, the q-values and the expected results.  For
the effective radius test and volume ratio tests, use the extended output
form, which checks each output of kernel.Fq. For 1-D tests, either specify
the q value or a list of q-values, and the corresponding I(q) value, or
list of I(q) values.

That is::

    tests = [
        [ {parameters}, q, I(q)],
        [ {parameters}, [q], [I(q)] ],
        [ {parameters}, [q1, q2, ...], [I(q1), I(q2), ...]],

        [ {parameters}, (qx, qy), I(qx, Iqy)],
        [ {parameters}, [(qx1, qy1), (qx2, qy2), ...],
                        [I(qx1, qy1), I(qx2, qy2), ...]],

        [ {parameters}, q, F(q), F^2(q), R_eff, V, V_r ],
        ...
    ]

Parameters are *key:value* pairs, where key is one of the parameters of the
model and value is the value to use for the test.  Any parameters not given
in the parameter list will take on the default parameter value.

Precision defaults to 5 digits (relative).
"""
from __future__ import print_function

import argparse
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

from .core import list_models, load_model_info, build_model
from .direct_model import call_kernel, call_Fq
from .exception import annotate_exception
from .modelinfo import expand_pars
from .kernelcl import use_opencl
from .kernelcuda import use_cuda
from . import product

# pylint: disable=unused-import
try:
    from typing import List, Iterator, Callable, Any, Dict, Tuple, Union
    from .modelinfo import ParameterTable, ParameterSet, TestCondition, ModelInfo
    from .kernel import KernelModel
    DType = Union[None, str, np.dtype]
except ImportError:
    pass
# pylint: enable=unused-import

def make_suite(loaders, models):
    # type: (List[str], List[str]) -> unittest.TestSuite
    """
    Construct the pyunit test suite.

    *loaders* is the list of kernel drivers to use (dll, opencl or cuda).
    For python model the python driver is always used.

    *models* is the list of models to test, or *["all"]* to test all models.
    """
    suite = unittest.TestSuite()

    try:
        # See if the first model parses as a model group
        group = list_models(models[0])
        skip = models[1:]
        models = group
    except Exception:
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
        test_name = "%s[python]"%model_info.name
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
            test_name = "%s[dll]"%model_info.name
            test_method_name = "test_%s_dll" % model_info.id
            test = ModelTestCase(test_name, model_info,
                                 test_method_name,
                                 platform="dll",
                                 dtype="double",
                                 stash=stash)
            suite.addTest(test)

        # test using opencl if desired and available
        if 'opencl' in loaders and use_opencl():
            test_name = "%s[opencl]"%model_info.name
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

        # test using cuda if desired and available
        if 'cuda' in loaders and use_cuda():
            test_name = "%s[cuda]" % model_info.id
            test_method_name = "test_%s_cuda" % model_info.id
            # Using dtype=None so that the models that are only
            # correct for double precision are not tested using
            # single precision.  The choice is determined by the
            # presence of *single=False* in the model file.
            test = ModelTestCase(test_name, model_info,
                                 test_method_name,
                                 platform="cuda", dtype=None,
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
            self._failures = []  # Set of failed target values

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
                # test that Fq will run, and return R_eff, V, V_r
                ({}, 0.1, None, None, None, None, None),
                ]
            tests = smoke_tests
            #tests = []
            if self.info.tests is not None:
                tests += self.info.tests
            S_tests = [test for test in tests if '@S' in test[0]]
            P_tests = [test for test in tests if '@S' not in test[0]]
            del self._failures[:]
            try:
                model = build_model(self.info, dtype=self.dtype,
                                    platform=self.platform)
                results = [self.run_one(model, test) for test in P_tests]
                for test in S_tests:
                    # pull the S model name out of the test defn
                    pars = test[0].copy()
                    s_name = pars.pop('@S')
                    ps_test = [pars] + list(test[1:])
                    #print("PS TEST PARAMS!!!",ps_test)
                    # build the P@S model
                    s_info = load_model_info(s_name)
                    ps_info = product.make_product_info(self.info, s_info)
                    ps_model = build_model(ps_info, dtype=self.dtype,
                                           platform=self.platform)
                    # run the tests
                    #self.info = ps_model.info
                    #print("SELF.INFO PARAMS!!!",[p.id for p in self.info.parameters.call_parameters])
                    #print("PS MODEL PARAMETERS:",[p.id for p in ps_model.info.parameters.call_parameters])
                    results.append(self.run_one(ps_model, ps_test))

                if self.stash:
                    for test, target, actual in zip(tests, self.stash[0], results):
                        assert np.all(abs(target-actual) < 5e-5*abs(actual)), \
                            ("GPU/CPU comparison expected %s but got %s for %s"
                             % (target, actual, test[0]))
                else:
                    self.stash.append(results)

                # Check for missing tests.  Only do so for the "dll" tests
                # to reduce noise from both opencl and cuda, and because
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

            if self._failures:
                msg = ("The following assertions failed for %s:\n  %s"
                       % (self.test_name, "\n  ".join(self._failures)))
                raise AssertionError(msg)

        def _find_missing_tests(self):
            # type: () -> None
            """make sure there are 1D and 2D tests as appropriate"""
            model_has_1D = True
            model_has_2D = any(p.type == 'orientation'
                               for p in self.info.parameters.kernel_parameters)

            # Lists of tests that have a result that is not None
            single = [test for test in self.info.tests
                      if not isinstance(test[2], list) and test[2] is not None]
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
            if model_has_1D and not (tests_has_1D_single or tests_has_1D_multiple):
                missing.append("1D")
            if model_has_2D and not (tests_has_2D_single or tests_has_2D_multiple):
                missing.append("2D")

            return missing

        def run_one(self, model, test):
            # type: (KernelModel, TestCondition) -> None
            """Run a single test case."""
            user_pars, x, y = test[:3]
            #print("PS MODEL PARAMETERS:",[p.id for p in model.info.parameters.call_parameters])
            pars = expand_pars(model.info.parameters, user_pars)
            invalid = invalid_pars(model.info.parameters, pars)
            if invalid:
                raise ValueError("Unknown parameters in test: " + ", ".join(invalid))

            if not isinstance(y, list):
                y = [y]
            if not isinstance(x, list):
                x = [x]

            self.assertEqual(len(y), len(x))

            if isinstance(x[0], tuple):
                qx, qy = zip(*x)
                q_vectors = [np.array(qx), np.array(qy)]
            else:
                q_vectors = [np.array(x)]

            kernel = model.make_kernel(q_vectors)
            if len(test) == 3 or len(test) == 4:
                actual = call_kernel(kernel, pars)
                self._check_vectors(x, y, actual, 'I')
                if len(test) == 4:
                    results = getattr(kernel, 'results', lambda: {})
                    self._check_struct(x, test[3], results())
                return actual
            elif len(test) == 7:
                y1 = y
                y2 = test[3] if isinstance(test[3], list) else [test[3]]
                F, Fsq, R_eff, volume, volume_ratio = call_Fq(kernel, pars)
                if F is not None:  # F is none for models with Iq instead of Fq
                    self._check_vectors(x, y1, F, 'F')
                self._check_vectors(x, y2, Fsq, 'F^2')
                self._check_scalar(test[4], R_eff, 'R_eff')
                self._check_scalar(test[5], volume, 'volume')
                self._check_scalar(test[6], volume_ratio, 'form:shell ratio')
                return Fsq
            else:
                self._failures.append('wrong number or results for %s => %s'
                                      % (str(user_pars), str(test[2:])))
                return None

        def _check_scalar(self, target, actual, name):
            if not is_near(target, actual, 5):
                self._failures.append('%s: expected:%s; actual:%s'
                                      % (name, target, actual))

        def _check_vectors(self, x, target, actual, name='I'):
            if not len(actual) > 0:
                self._failures.append('%s(...) expected return value'%name)
                return
            if target is None:
                return
            if len(target) != len(actual):
                self._failures.append('%s(...) returned wrong length in %s'
                                      % (name, str(actual)))
                return
            for xi, yi, actual_yi in zip(x, target, actual):
                if not is_near(yi, actual_yi, 5):
                    # convert array to list so we have comma separated output
                    actual = actual.tolist()
                    self._failures.append('%s(%s): expected:%s; actual:%s'
                                          % (name, xi, target, actual))
                    # Whole list is printed on any error, so fail on first
                    return

        def _check_struct(self, x, target, actual):
            for k, v in target.items():
                if k not in actual:
                    self._failures.append('key %r not in returned value' % k)
                    continue
                target_k = v
                actual_k = actual[k]
                if np.isscalar(actual_k):  # number
                    self._check_scalar(target_k, actual_k, k)
                elif isinstance(actual_k, np.ndarray): # vector
                    self._check_vectors(x, target_k, actual_k, k)
                elif isinstance(actual_k, tuple) and len(actual_k) == 2:
                    # Intermediate is returned as (Q, I(Q)) pair.  The test
                    # is set up to use |Q| or (Qx, Qy) as input, but the
                    # returned result has (|Q|,) or (Qx, Qy).  For now we will
                    # ignore the first result.
                    #self._check_vectors(x, x, actual_k[0], k+" (Q)")
                    self._check_vectors(x, target_k, actual_k[1], k)
                elif isinstance(actual_k, dict):
                    self._check_struct(x, target_k, actual_k)
                else:
                    self._failures.append('key %s: unexpected value %r'
                                          % (k, actual_k))
            # Ignore effective_radius in returned result; it will be dropped
            actual.pop('effective_radius', None)
            if any(k not in target for k in actual):
                self._failures.append('intermediate results missing:')
                for k, v in actual.items():
                    if k not in target:
                        if isinstance(v, tuple) and len(v) == 2:
                            v = v[1].tolist()
                        elif isinstance(v, np.ndarray):
                            v = v.tolist()
                        self._failures.append('  "%s": %r,' % (k, v))

    return ModelTestCase

def invalid_pars(partable, pars):
    # type: (ParameterTable, Dict[str, float]) -> List[str]
    """
    Return a list of parameter names that are not part of the model.
    """
    names = set(p.id for p in partable.call_parameters)
    invalid = []
    for par in sorted(pars.keys()):
        # Ignore the R_eff mode parameter when checking for valid parameters.
        # It is an allowed parameter for a model even though it does not exist
        # in the parameter table.  The call_Fq() function pops it from the
        # parameter list and sends it directly to kernel.Fq().
        if par == product.RADIUS_MODE_ID:
            continue
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

    *taget* zero and inf should match *actual* zero and inf.  If you want to
    accept eps for zero, choose a value such as 1e-10, which must match up to
    +/- 1e-15 when *digits* is the default value of 5.

    If *target* is None, then just make sure that *actual* is not NaN.

    If *target* is NaN, make sure *actual* is NaN.
    """
    if target is None:
        # target is None => actual cannot be NaN
        return not np.isnan(actual)
    elif target == 0.:
        # target is 0. => actual must be 0.
        # Note: if small values are allowed, then use maybe test zero against eps instead?
        return actual == 0.
    elif np.isfinite(target):
        shift = np.ceil(np.log10(abs(target)))
        return abs(target-actual) < 1.5*10**(shift-digits)
    elif target == actual:
        # target is inf => actual must be inf of same sign
        return True
    else:
        # target is NaN => actual must be NaN
        return np.isnan(target) == np.isnan(actual)

# CRUFT: old interface; should be deprecated and removed
def run_one(model_name):
    # type: (str) -> str
    """
    [Deprecated] Run the tests associated with *model_name*.

    Use the following instead::

        succss, output = check_model(load_model_info(model_name))
    """
    # msg = "use check_model(model_info) rather than run_one(model_name)"
    # warnings.warn(msg, category=DeprecationWarning, stacklevel=2)
    try:
        model_info = load_model_info(model_name)
    except Exception:
        output = traceback.format_exc()
        return output

    _, output = check_model(model_info)
    return output

def check_model(model_info):
    # type: (ModelInfo) -> Tuple[bool, str]
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
    loaders = ['opencl' if use_opencl() else 'cuda' if use_cuda() else 'dll']
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


def model_tests():
    # type: () -> Iterator[Callable[[], None]]
    """
    Test runner visible to nosetests.

    Run "nosetests sasmodels" on the command line to invoke it.
    """
    loaders = ['dll']
    if use_opencl():
        loaders.append('opencl')
    if use_cuda():
        loaders.append('cuda')
    tests = make_suite(loaders, ['all'])
    def _build_test(test):
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
        yield _build_test(test)


def main():
    # type: () -> int
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

    parser = argparse.ArgumentParser(description="Test SasModels Models")
    parser.add_argument("-v", "--verbose", action="store_const",
                        default=1, const=2, help="Use verbose output")
    parser.add_argument("-e", "--engine", default="all",
                        help="Engines on which to run the test.  "
                        "Valid values are opencl, cuda, dll, and all. "
                        "Defaults to all if no value is given")
    parser.add_argument("models", nargs="*",
                        help="The names of the models to be tested.  "
                        "If the first model is 'all', then all but the listed "
                        "models will be tested.  See core.list_models() for "
                        "names of other groups, such as 'py' or 'single'.")
    opts = parser.parse_args()

    if opts.engine == "opencl":
        if not use_opencl():
            print("opencl is not available")
            return 1
        loaders = ['opencl']
    elif opts.engine == "dll":
        loaders = ["dll"]
    elif opts.engine == "cuda":
        if not use_cuda():
            print("cuda is not available")
            return 1
        loaders = ['cuda']
    elif opts.engine == "all":
        loaders = ['dll']
        if use_opencl():
            loaders.append('opencl')
        if use_cuda():
            loaders.append('cuda')
    else:
        print("unknown engine " + opts.engine)
        return 1

    runner = TestRunner(verbosity=opts.verbose, **test_args)
    result = runner.run(make_suite(loaders, opts.models))
    return 1 if result.failures or result.errors else 0


if __name__ == "__main__":
    sys.exit(main())
