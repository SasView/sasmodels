"""
GPU data alignment.

Some web sites say that maximizing performance for OpenCL code requires
aligning data on certain memory boundaries.

:func:`data_alignment` queries all devices in the OpenCL context, returning
the most restrictive alignment.

:func:`align_data` aligns an existing array.

:func:`align_empty` to create a new array of the correct alignment.

Note:  This code is unused. So far, tests have not demonstrated any
improvement from forcing correct alignment.  The tests should
be repeated with arrays forced away from the target boundaries
to decide whether it is really required.
"""

import numpy as np
import pyopencl as cl

def data_alignment(context):
    """
    Return the desired byte alignment across all devices.
    """
    # Note: rely on the fact that alignment will be 2^k
    return max(d.min_data_type_align_size for d in context.devices)

def align_empty(shape, dtype, alignment=128):
    size = np.prod(shape)
    dtype = np.dtype(dtype)
    # allocate array with extra space for alignment
    extra = alignment//dtype.itemsize - 1
    result = np.empty(size+extra, dtype)
    # build a view into allocated array which starts on a boundary
    offset = (result.ctypes.data%alignment)//dtype.itemsize
    view = np.reshape(result[offset:offset+size], shape)
    return view

def align_data(v, dtype, alignment=128):
    # if v is contiguous, aligned, and of the correct type then just return v
    view = align_empty(v.shape, dtype, alignment=alignment)
    view[:] = v
    return view

def work_group_size(context, kernel):
    """
    Return the desired work group size for the context.
    """
    max(kernel.get_work_group_info(cl.kernel_work_group_info.PREFERRED_WORK_GROUP_SIZE_MULTIPLE, d)
        for d in context.devices)


