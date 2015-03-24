import numpy as np
from numpy import pi, cos

from .generate import F64

class PyModel(object):
    def __init__(self, info):
        self.info = info

    def __call__(self, q_input):
        kernel = self.info['Iqxy'] if q_input.is_2D else self.info['Iq']
        return PyKernel(kernel, self.info, q_input)

    # pylint: disable=no-self-use
    def make_input(self, q_vectors):
        return PyInput(q_vectors, dtype=F64)

    def release(self):
        pass

class PyInput(object):
    """
    Make q data available to the gpu.

    *q_vectors* is a list of q vectors, which will be *[q]* for 1-D data,
    and *[qx, qy]* for 2-D data.  Internally, the vectors will be reallocated
    to get the best performance on OpenCL, which may involve shifting and
    stretching the array to better match the memory architecture.  Additional
    points will be evaluated with *q=1e-3*.

    *dtype* is the data type for the q vectors. The data type should be
    set to match that of the kernel, which is an attribute of
    :class:`GpuProgram`.  Note that not all kernels support double
    precision, so even if the program was created for double precision,
    the *GpuProgram.dtype* may be single precision.

    Call :meth:`release` when complete.  Even if not called directly, the
    buffer will be released when the data object is freed.
    """
    def __init__(self, q_vectors, dtype):
        self.nq = q_vectors[0].size
        self.dtype = dtype
        self.is_2D = (len(q_vectors) == 2)
        self.q_vectors = [np.ascontiguousarray(q, self.dtype) for q in q_vectors]
        self.q_pointers = [q.ctypes.data for q in self.q_vectors]

    def release(self):
        self.q_vectors = []

class PyKernel(object):
    """
    Callable SAS kernel.

    *kernel* is the DllKernel object to call.

    *info* is the module information

    *q_input* is the DllInput q vectors at which the kernel should be
    evaluated.

    The resulting call method takes the *pars*, a list of values for
    the fixed parameters to the kernel, and *pd_pars*, a list of (value,weight)
    vectors for the polydisperse parameters.  *cutoff* determines the
    integration limits: any points with combined weight less than *cutoff*
    will not be calculated.

    Call :meth:`release` when done with the kernel instance.
    """
    def __init__(self, kernel, info, q_input):
        self.info = info
        self.q_input = q_input
        self.res = np.empty(q_input.nq, q_input.dtype)
        dim = '2d' if q_input.is_2D else '1d'
        # Loop over q unless user promises that the kernel is vectorized by
        # taggining it with vectorized=True
        if not getattr(kernel, 'vectorized', False):
            if dim == '2d':
                def vector_kernel(qx, qy, *args):
                    return np.array([kernel(qxi, qyi, *args)
                                     for qxi, qyi in zip(qx, qy)])
            else:
                def vector_kernel(q, *args):
                    return np.array([kernel(qi, *args)
                                     for qi in q])
            self.kernel = vector_kernel
        else:
            self.kernel = kernel
        fixed_pars = info['partype']['fixed-' + dim]
        pd_pars = info['partype']['pd-' + dim]
        vol_pars = info['partype']['volume']

        # First two fixed pars are scale and background
        pars = [p[0] for p in info['parameters'][2:]]
        offset = len(self.q_input.q_vectors)
        self.args = self.q_input.q_vectors + [None] * len(pars)
        self.fixed_index = np.array([pars.index(p) + offset
                                     for p in fixed_pars[2:]])
        self.pd_index = np.array([pars.index(p) + offset
                                  for p in pd_pars])
        self.vol_index = np.array([pars.index(p) + offset
                                   for p in vol_pars])
        try: self.theta_index = pars.index('theta') + offset
        except ValueError: self.theta_index = -1

        # Caller needs fixed_pars and pd_pars
        self.fixed_pars = fixed_pars
        self.pd_pars = pd_pars

    def __call__(self, fixed, pd, cutoff=1e-5):
        #print "fixed",fixed
        #print "pd", pd
        args = self.args[:]  # grab a copy of the args
        form, form_volume = self.kernel, self.info['form_volume']
        # First two fixed
        scale, background = fixed[:2]
        for index, value in zip(self.fixed_index, fixed[2:]):
            args[index] = float(value)
        res = _loops(form, form_volume, cutoff, scale, background, args,
                     pd, self.pd_index, self.vol_index, self.theta_index)

        return res

    def release(self):
        self.q_input = None

def _loops(form, form_volume, cutoff, scale, background,
           args, pd, pd_index, vol_index, theta_index):
    """
    *form* is the name of the form function, which should be vectorized over
    q, but otherwise have an interface like the opencl kernels, with the
    q parameters first followed by the individual parameters in the order
    given in model.parameters (see :mod:`sasmodels.generate`).

    *form_volume* calculates the volume of the shape.  *vol_index* gives
    the list of volume parameters

    *cutoff* ignores the corners of the dispersion hypercube

    *scale*, *background* multiplies the resulting form and adds an offset

    *args* is the prepopulated set of arguments to the form function, starting
    with the q vectors, and including slots for all the parameters.  The
    values for the parameters will be substituted with values from the
    dispersion functions.

    *pd* is the list of dispersion parameters

    *pd_index* are the indices of the dispersion parameters in *args*

    *vol_index* are the indices of the volume parameters in *args*

    *theta_index* is the index of the theta parameter for the sasview
    spherical correction, or -1 if there is no angular dispersion
    """

    ################################################################
    #                                                              #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #   !!                                                    !!   #
    #   !!  KEEP THIS CODE CONSISTENT WITH KERNEL_TEMPLATE.C  !!   #
    #   !!                                                    !!   #
    #   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   #
    #                                                              #
    ################################################################

    weight = np.empty(len(pd), 'd')
    if weight.size > 0:
        # weight vector, to be populated by polydispersity loops

        # identify which pd parameters are volume parameters
        vol_weight_index = np.array([(index in vol_index) for index in pd_index])

        # Sort parameters in decreasing order of pd length
        Npd = np.array([len(pdi[0]) for pdi in pd], 'i')
        order = np.argsort(Npd)[::-1]
        stride = np.cumprod(Npd[order])
        pd = [pd[index] for index in order]
        pd_index = pd_index[order]
        vol_weight_index = vol_weight_index[order]

        fast_value = pd[0][0]
        fast_weight = pd[0][1]
    else:
        stride = np.array([1])
        vol_weight_index = slice(None, None)
        # keep lint happy
        fast_value = [None]
        fast_weight = [None]

    ret = np.zeros_like(args[0])
    norm = np.zeros_like(ret)
    vol = np.zeros_like(ret)
    vol_norm = np.zeros_like(ret)
    for k in range(stride[-1]):
        # update polydispersity parameter values
        fast_index = k % stride[0]
        if fast_index == 0:  # bottom loop complete ... check all other loops
            if weight.size > 0:
                for i, index, in enumerate(k % stride):
                    args[pd_index[i]] = pd[i][0][index]
                    weight[i] = pd[i][1][index]
        else:
            args[pd_index[0]] = fast_value[fast_index]
            weight[0] = fast_weight[fast_index]
        # This computes the weight, and if it is sufficient, calls the
        # scattering function and adds it to the total.  If there is a volume
        # normalization, it will also be added here.
        # Note: make sure this is consistent with the code in PY_LOOP_BODY!!
        # Note: can precompute w1*w2*...*wn
        w = np.prod(weight)
        if w > cutoff:
            I = form(*args)
            positive = (I >= 0.0)

            # Note: can precompute spherical correction if theta_index is not
            # the fast index. Correction factor for spherical integration
            #spherical_correction = abs(cos(pi*args[phi_index])) if phi_index>=0 else 1.0
            spherical_correction = (abs(cos(pi * args[theta_index])) * pi / 2
                                    if theta_index >= 0 else 1.0)
            #spherical_correction = 1.0
            ret += w * I * spherical_correction * positive
            norm += w * positive

            # Volume normalization.
            # If there are "volume" polydispersity parameters, then these
            # will be used to call the form_volume function from the user
            # supplied kernel, and accumulate a normalized weight.
            # Note: can precompute volume norm if fast index is not a volume
            if form_volume:
                vol_args = [args[index] for index in vol_index]
                vol_weight = np.prod(weight[vol_weight_index])
                vol += vol_weight * form_volume(*vol_args) * positive
                vol_norm += vol_weight * positive

    positive = (vol * vol_norm != 0.0)
    ret[positive] *= vol_norm[positive] / vol[positive]
    result = scale * ret / norm + background
    return result
