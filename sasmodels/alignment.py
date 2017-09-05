"""
GPU data alignment.

Some web sites say that maximizing performance for OpenCL code requires
aligning data on certain memory boundaries.  The following functions
provide this service:

:func:`align_data` aligns an existing array, returning a new array of the
correct alignment.

:func:`align_empty` to create an empty array of the correct alignment.

Set alignment to :func:`gpu.environment()` attribute *boundary*.

Note:  This code is unused. So far, tests have not demonstrated any
improvement from forcing correct alignment.  The tests should
be repeated with arrays forced away from the target boundaries
to decide whether it is really required.
"""
import numpy as np  # type: ignore

def align_empty(shape, dtype, alignment=128):
    """
    Return an empty array aligned on the alignment boundary.
    """
    size = np.prod(shape)
    dtype = np.dtype(dtype)
    # allocate array with extra space for alignment
    extra = alignment//dtype.itemsize - 1
    result = np.empty(size+extra, dtype)
    # build a view into allocated array which starts on a boundary
    offset = (result.ctypes.data%alignment)//dtype.itemsize
    view = np.reshape(result[offset:offset+size], shape)
    return view

def align_data(x, dtype, alignment=128):
    """
    Return a copy of an array on the alignment boundary.
    """
    # if x is contiguous, aligned, and of the correct type then just return x
    view = align_empty(x.shape, dtype, alignment=alignment)
    view[:] = x
    return view
