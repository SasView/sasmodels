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
from os.path import basename, splitext, join as joinpath, exists, dirname

try:
    # Python 3.5 and up
    from importlib.util import spec_from_file_location, module_from_spec  # type: ignore
    def load_module_from_path(fullname, path):
        # type: (str, str) -> "module"
        """load module from *path* as *fullname*"""
        spec = spec_from_file_location(fullname, os.path.expanduser(path))
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
except ImportError:
    # CRUFT: python 2
    import imp
    def load_module_from_path(fullname, path):
        # type: (str, str) -> "module"
        """load module from *path* as *fullname*"""
        # Clear out old definitions, if any
        if fullname in sys.modules:
            del sys.modules[fullname]
        if path.endswith(".py") and os.path.exists(path) and os.path.exists(path+"c"):
            # remove automatic pyc file before loading a py file
            os.unlink(path+"c")
        module = imp.load_source(fullname, os.path.expanduser(path))
        return module

_MODULE_CACHE = {} # type: Dict[str, Tuple("module", int)]
_MODULE_DEPENDS = {} # type: Dict[str, List[str]]
_MODULE_DEPENDS_STACK = [] # type: List[str]
def load_custom_kernel_module(path):
    # type: str -> "module"
    """load SAS kernel from *path* as *sasmodels.custom.modelname*"""
    # Pull off the last .ext if it exists; there may be others
    name = basename(splitext(path)[0])
    path = os.path.expanduser(path)

    # Reload module if necessary.
    if need_reload(path):
        # Assume the module file is the only dependency
        _MODULE_DEPENDS[path] = set([path])

        # Load the module while pushing it onto the dependency stack.  If
        # this triggers any submodules, then they will add their dependencies
        # to this module as the "working_on" parent.  Pop the stack when the
        # module is loaded.
        _MODULE_DEPENDS_STACK.append(path)
        module = load_module_from_path('sasmodels.custom.'+name, path)
        _MODULE_DEPENDS_STACK.pop()

        # Include external C code in the dependencies.  We are looking
        # for module.source and assuming that it is a list of C source files
        # relative to the module itself.  Any files that do not exist,
        # such as those in the standard libraries, will be ignored.
        # TODO: look in builtin module path for standard c sources
        # TODO: share code with generate.model_sources
        c_sources = getattr(module, 'source', None)
        if isinstance(c_sources, (list, tuple)):
            _MODULE_DEPENDS[path].update(_find_sources(path, c_sources))

        # Cache the module, and tag it with the newest timestamp
        timestamp = max(os.path.getmtime(f) for f in _MODULE_DEPENDS[path])
        _MODULE_CACHE[path] = module, timestamp

        #print("loading", os.path.basename(path), _MODULE_CACHE[path][1],
        #    [os.path.basename(p) for p in _MODULE_DEPENDS[path]])

    # Add path and all its dependence to the parent module, if there is one.
    if _MODULE_DEPENDS_STACK:
        working_on = _MODULE_DEPENDS_STACK[-1]
        _MODULE_DEPENDS[working_on].update(_MODULE_DEPENDS[path])

    return _MODULE_CACHE[path][0]

def need_reload(path):
    # type: str -> bool
    """
    Return True if any path dependencies have a timestamp newer than the time
    when the path was most recently loaded.
    """
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
    #print("reload", any(cache_time < os.path.getmtime(p) for p in depends))
    #for f in depends: print(">>>  ", f, os.path.getmtime(f))
    return any(cache_time < os.path.getmtime(p) for p in depends)

def _find_sources(path, source_list):
    # type: (str, List[str]) -> List[str]
    """
    Return a list of the sources relative to base file; ignore any that
    are not found.
    """
    root = dirname(path)
    found = []
    for source_name in source_list:
        source_path = joinpath(root, source_name)
        if exists(source_path):
            found.append(source_path)
    return found
