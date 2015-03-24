from astroid import MANAGER
from astroid import nodes

def register(linter):
    MANAGER.register_transform(nodes.Module, transform)

def transform(module):
    #print("processing",module.name)
    if module.name == 'pyopencl':
        import pyopencl
