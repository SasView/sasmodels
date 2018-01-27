"""
py.test hooks for sasmodels

Hooks for running sasmodels tests via py.test.

*pytest_collection_modifyitems* adds the test description to the end of
the test name.  This is needed for the generated list of tests is sasmodels,
where each test has a description giving the name of the model.  For example
"model_tests::[3]" becomes "model_tests::[3]::bcc_paracrystal-dll".  Need to
leave the "::[3]" in the name since that is the way you can indicate this
test specifically from the py.test command line. [This is perhaps because
the modifyitems hook is only called after test selection.]
"""
from __future__ import print_function

import pytest
from _pytest.unittest import TestCaseFunction

USE_DOCSTRING_AS_DESCRIPTION = True
def pytest_collection_modifyitems(session, config, items):
    """
    Add description to the test node id if item is a function and function
    has a description attribute or __doc__ attribute.
    """
    for item in items:
        #print(item.nodeid, type(item))
        #for attr in dir(item): not attr.startswith('__') and print(attr, getattr(item, attr))
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
                description = str(item._args) if len(item._args) > 1 else str(item._args[0])
            #print(item.nodeid, description, item._args)

            if description:
                # Strip spaces from start and end and strip dots from end
                # pytest converts '.' to '::' on output for some reason.
                description = description.strip().rstrip('.')
                # Join multi-line descriptions into a single line
                if '\n' in description:
                    description = " ".join(line.strip() for line in description.split('\n'))

            # Set the description as part of the node identifier.
            if description:
                #print(type(item), dir(item))
                #print(item.nodeid, description)
                #print(item.location, item.name, item.nodeid, item.originalname)
                # Note: leave the current name mostly as-is since the prefix
                # is needed to specify the nth test from a list of tests.
                #print("updating with", description)
                item.name += "::" + description
            #print("=>", item.nodeid)
