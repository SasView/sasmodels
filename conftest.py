"""
pytest hooks for sasmodels

Hooks for running sasmodels tests via pytest.

*pytest_collection_modifyitems* adds the test description to the end of
the test name.  This is needed for the generated list of tests is sasmodels,
where each test has a description giving the name of the model.  For example
"model_tests::[3]" becomes "model_tests::[3]::bcc_paracrystal-dll".  Need to
leave the "::[3]" in the name since that is the way you can indicate this
test specifically from the pytest command line. [This is perhaps because
the modifyitems hook is only called after test selection.]

*pytest_ignore_collect* skips kernelcl.py if pyopencl cannot be imported.
"""
from __future__ import print_function

import os.path
import inspect

import pytest
from _pytest.unittest import TestCaseFunction

# CRUFT: pytest<5.4 does not have from_parent class method.
function_test = getattr(pytest.Function, 'from_parent', pytest.Function)

def pytest_pycollect_makeitem(collector, name, obj):
    """
    Convert test generator into list of function tests so that pytest doesn't
    complain about deprecated yield tests.

    Note that unlike nose, the tests are generated and saved instead of run
    immediately.  This means that any dynamic context, such as a for-loop
    variable, must be captured by wrapping the yield result in a function call.

    For example::

        for value in 1, 2, 3:
            for test in test_cases:
                yield test, value

    will need to be changed to::

        def build_test(test, value):
            return test, value
        for value in 1, 2, 3:
            for test in test_cases:
                yield build_test(test, value)

    This allows the context (test and value) to be captured by lexical closure
    in build_test. See https://stackoverflow.com/a/233835/6195051.
    """
    if collector.istestfunction(obj, name) and is_generator(obj):
        tests = []
        for number, yielded in enumerate(obj()):
            index, call = split_yielded_test(yielded, number)
            description = getattr(call, 'description', name+index)
            test = function_test(
                parent=collector, name=description, callobj=call)
            tests.append(test)
        return tests

def is_generator(func):
    """
    Returns True if function has yield.
    """
    # Cribbed from _pytest.compat is_generator and iscoroutinefunction; these
    # may not be available as part of pytest 4.
    coroutine = (getattr(func, '_is_coroutine', False) or
                 getattr(inspect, 'iscoroutinefunction', lambda f: False)(func))
    generator = inspect.isgeneratorfunction(func)
    return generator and not coroutine

def split_yielded_test(obj, number):
    if not isinstance(obj, (tuple, list)):
        obj = (obj,)
    assert callable(obj[0])
    index = "[%d]"%number
    call = obj[0]
    return index, call

USE_DOCSTRING_AS_DESCRIPTION = True
def pytest_collection_modifyitems(session, config, items):
    """
    Add description to the test node id if item is a function and function
    has a description attribute or __doc__ attribute.
    """
    for item in items:
        if isinstance(item, pytest.Function):
            if isinstance(item, TestCaseFunction):
                # TestCase uses item.name to find the method so skip
                continue
            function = item.obj

            # If the test case provides a "description" attribute then use it
            # as an extended description.  If there is no description attribute,
            # then perhaps use the test docstring.
            if USE_DOCSTRING_AS_DESCRIPTION:
                description = getattr(function, 'description', function.__doc__)
            else:
                description = getattr(function, 'description', "")

            # If description is not supplied but yield args are, then use the
            # yield args for the description
            if not description and getattr(item, '_args', ()):
                description = (str(item._args) if len(item._args) > 1
                               else str(item._args[0]))

            # Set the description as part of the node identifier.
            if description:
                # Strip spaces from start and end and strip dots from end
                # pytest converts '.' to '::' on output for some reason.
                description = description.strip().rstrip('.')
                # Join multi-line descriptions into a single line
                if '\n' in description:
                    description = " ".join(line.strip()
                                           for line in description.split('\n'))

                # Note: leave the current name mostly as-is since the prefix
                # is needed to specify the nth test from a list of tests.
                item.name += "::" + description
