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

_MODULE_CACHE = {}
_MODULE_DEPENDS = {}
_MODULE_DEPENDS_STACK = []
def load_custom_kernel_module(path):
    """load SAS kernel from *path* as *sasmodels.custom.modelname*"""
    # Pull off the last .ext if it exists; there may be others
    name = basename(splitext(path)[0])
    path = os.path.expanduser(path)

    # reload module if necessary
    if need_reload(path):
        # Push to the next dependency level
        _MODULE_DEPENDS_STACK.append(path)
        _MODULE_DEPENDS[path] = set([path])

        # Load module into the 'sasmodels.custom' name space.
        # If this triggers any submodule loads then they will be added
        # as dependencies below when _MODULE_DEPENDS_STACK is not empty.
        module = load_module_from_path('sasmodels.custom.'+name, path)

        # Pop the dependency level
        _MODULE_DEPENDS_STACK.pop()

        # TODO: include external C code in the dependencies
        # If we had the model info structure we could do the following:
        #    _MODEL_DEPENDS[path].extend(generate.model_sources(info))
        # but at this point all we have is the module.  Don't want to
        # repeat the logic in modelinfo.make_model_info.

        # Cache the module with the newest timestamp
        timestamp = max(os.path.getmtime(f) for f in _MODULE_DEPENDS[path])
        _MODULE_CACHE[path] = module, timestamp

        #print("loading", os.path.basename(path), _MODULE_CACHE[path][1],
        #    [os.path.basename(p) for p in _MODULE_DEPENDS[path]])

    if _MODULE_DEPENDS_STACK:
        # Add child and all its dependence to the parent module
        working_on = _MODULE_DEPENDS_STACK[-1]
        _MODULE_DEPENDS[working_on].update(_MODULE_DEPENDS[path])

    return _MODULE_CACHE[path][0]

def need_reload(path):
    # TODO: fails if a dependency has a modification time in the future
    # If the newest dependency has a time stamp in the future, then this
    # will be recorded as the cached time.  When a second dependency
    # is updated to the current time stamp, it will still be considered
    # older than the current build and the reload will not be triggered.
    # Could instead treat all future times as 0 here and in the code above
    # which records the newest timestamp.  This will force a reload when
    # the future time is reached, but other than that should perform
    # correctly.  Probably not worth the extra code...
    _, cache_time = _MODULE_CACHE.get(path, (None, -1))
    depends = _MODULE_DEPENDS.get(path, [path])
    return any(cache_time < os.path.getmtime(p) for p in depends)
