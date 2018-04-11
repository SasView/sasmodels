"""
Custom Models
-------------

This is a place holder for the custom models namespace.  When models are
loaded from a file by :func:`generate.load_kernel_module` they are loaded
as if they exist in *sasmodels.custom*.  This package needs to exist for this
to occur without error.
"""
from __future__ import division, print_function

import sys
import os
from os.path import basename, splitext

try:
    # Python 3.5 and up
    from importlib.util import spec_from_file_location, module_from_spec  # type: ignore
    def load_module_from_path(fullname, path):
        """load module from *path* as *fullname*"""
        spec = spec_from_file_location(fullname, os.path.expanduser(path))
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
except ImportError:
    # CRUFT: python 2
    import imp
    def load_module_from_path(fullname, path):
        """load module from *path* as *fullname*"""
        # Clear out old definitions, if any
        if fullname in sys.modules:
            del sys.modules[fullname]
        if path.endswith(".py") and os.path.exists(path) and os.path.exists(path+"c"):
            # remove automatic pyc file before loading a py file
            os.unlink(path+"c")
        module = imp.load_source(fullname, os.path.expanduser(path))
        return module

def load_custom_kernel_module(path):
    """load SAS kernel from *path* as *sasmodels.custom.modelname*"""
    # Pull off the last .ext if it exists; there may be others
    name = basename(splitext(path)[0])
    # Placing the model in the 'sasmodels.custom' name space.
    kernel_module = load_module_from_path('sasmodels.custom.'+name,
                                          os.path.expanduser(path))
    return kernel_module
