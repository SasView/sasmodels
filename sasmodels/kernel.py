"""
Execution kernel interface
==========================

:class:`KernelModel` defines the interface to all kernel models.
In particular, each model should provide a :meth:`KernelModel.make_kernel`
call which returns an executable kernel, :class:`Kernel`, that operates
on the given set of *q_vector* inputs.  On completion of the computation,
the kernel should be released, which also releases the inputs.
"""

from __future__ import division, print_function

# pylint: disable=unused-import
try:
    from typing import List, Any
except ImportError:
    pass
else:
    import numpy as np
    from .details import CallDetails
    from .modelinfo import ModelInfo
# pylint: enable=unused-import


class KernelModel(object):
    """
    Model definition for the compute engine.
    """
    info = None  # type: ModelInfo
    dtype = None # type: np.dtype
    def make_kernel(self, q_vectors):
        # type: (List[np.ndarray]) -> "Kernel"
        """
        Instantiate a kernel for evaluating the model at *q_vectors*.
        """
        raise NotImplementedError("need to implement make_kernel")

    def release(self):
        # type: () -> None
        """
        Free resources associated with the kernel.
        """
        #print("null release model")
        pass


class Kernel(object):
    """
    Instantiated model for the compute engine, applied to a particular *q*.

    Subclasses should define *__init__()* to set up the kernel inputs, and
    *_call_kernel()* to evaluate the kernel::

        def __init__(self, ...):
            ...
            self.q_input = <q-value class with nq attribute>
            self.info = <ModelInfo object>
            self.dim = <'1d' or '2d'>
            self.dtype = <kernel.dtype>
            size = 2*self.q_input.nq+4 if self.info.have_Fq else self.q_input.nq+4
            size = size + <extra padding if needed for kernel>
            self.result = np.empty(size, dtype=self.dtype)

        def _call_kernel(self, call_details, values, cutoff, magnetic,
                        radius_effective_mode):
            # type: (CallDetails, np.ndarray, np.ndarray, float, bool, int) -> None
            ... # call <kernel>
            nq = self.q_input.nq
            if self.info.have_Fq:  # models that compute both F and F^2
                end = 2*nq if have_Fq else nq
                self.result[0:end:2] = F**2
                self.result[1:end:2] = F
            else:
                end = nq
                self.result[0:end] = Fsq
            self.result[end + 0] = total_weight
            self.result[end + 1] = form_volume
            self.result[end + 2] = shell_volume
            self.result[end + 3] = radius_effective
    """
    #: Kernel dimension, either "1d" or "2d".
    dim = None  # type: str
    #: Model info.
    info = None  # type: ModelInfo
    #: Numerical precision for the computation.
    dtype = None  # type: np.dtype
    #: Q values at which the kernel is to be evaluated.
    q_input = None  # type: Any
    #: Place to hold result of *_call_kernel()* for subclass.
    result = None # type: np.ndarray

    def Iq(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool) -> np.ndarray
        r"""
        Returns I(q) from the polydisperse average scattering.

        .. math::

            I(q) = \text{scale} \cdot P(q) + \text{background}

        With the correct choice of model and contrast, setting *scale* to
        the volume fraction $V_f$ of particles should match the measured
        absolute scattering.  Some models (e.g., vesicle) have volume fraction
        built into the model, and do not need an additional scale.
        """
        _, F2, _, shell_volume, _ = self.Fq(call_details, values, cutoff,
                                            magnetic, radius_effective_mode=0)
        combined_scale = values[0]/shell_volume
        background = values[1]
        return combined_scale*F2 + background
    __call__ = Iq

    def Fq(self, call_details, values, cutoff, magnetic,
           radius_effective_mode=0):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool, int) -> np.ndarray
        r"""
        Returns <F(q)>, <F(q)^2>, effective radius, shell volume and
        form:shell volume ratio.  The <F(q)> term may be None if the
        form factor does not support direct computation of $F(q)$

        $P(q) = <F^2(q)>/<V>$ is used for structure factor calculations,

        .. math::

            I(q) = \text{scale} \cdot P(q) \cdot S(q) + \text{background}

        For the beta approximation, this becomes

        .. math::

            I(q) = \text{scale} P (1 + <F>^2/<F^2> (S - 1)) + \text{background}
                 = \text{scale}/<V> (<F^2> + <F>^2 (S - 1)) + \text{background}

        $<F(q)>$ and $<F^2(q)>$ are averaged by polydispersity in shape
        and orientation, with each configuration $x_k$ having form factor
        $F(q, x_k)$, weight $w_k$ and volume $V_k$.  The result is:

        .. math::

            P(q)=\frac{\sum w_k F^2(q, x_k) / \sum w_k}{\sum w_k V_k / \sum w_k}

        The form factor itself is scaled by volume and contrast to compute the
        total scattering.  This is then squared, and the volume weighted
        F^2 is then normalized by volume F.  For a given density, the number
        of scattering centers is assumed to scale linearly with volume.  Later
        scaling the resulting $P(q)$ by the volume fraction of particles
        gives the total scattering on an absolute scale. Most models
        incorporate the volume fraction into the overall scale parameter.  An
        exception is vesicle, which includes the volume fraction parameter in
        the model itself, scaling $F$ by $\surd V_f$ so that the math for
        the beta approximation works out.

        By scaling $P(q)$ by total weight $\sum w_k$, there is no need to make
        sure that the polydisperisity distributions normalize to one.  In
        particular, any distibution values $x_k$ outside the valid domain
        of $F$ will not be included, and the distribution will be implicitly
        truncated.  This is controlled by the parameter limits defined in the
        model (which truncate the distribution before calling the kernel) as
        well as any region excluded using the *INVALID* macro defined within
        the model itself.

        The volume used in the polydispersity calculation is the form volume
        for solid objects or the shell volume for hollow objects.  Shell
        volume should be used within $F$ so that the normalizing scale
        represents the volume fraction of the shell rather than the entire
        form.  This corresponds to the volume fraction of shell-forming
        material added to the solvent.

        The calculation of $S$ requires the effective radius and the
        volume fraction of the particles.  The model can have several
        different ways to compute effective radius, with the
        *radius_effective_mode* parameter used to select amongst them.  The
        volume fraction of particles should be determined from the total
        volume fraction of the form, not just the shell volume fraction.
        This makes a difference for hollow shapes, which need to scale
        the volume fraction by the returned volume ratio when computing $S$.
        For solid objects, the shell volume is set to the form volume so
        this scale factor evaluates to one and so can be used for both
        hollow and solid shapes.
        """
        self._call_kernel(call_details, values, cutoff, magnetic,
                          radius_effective_mode)
        #print("returned",self.q_input.q, self.result)
        nout = 2 if self.info.have_Fq and self.dim == '1d' else 1
        total_weight = self.result[nout*self.q_input.nq + 0]
        # Note: total_weight = sum(weight > cutoff), with cutoff >= 0, so it
        # is okay to test directly against zero.  If weight is zero then I(q),
        # etc. must also be zero.
        if total_weight == 0.:
            total_weight = 1.
        # Note: shell_volume == form_volume for solid objects
        form_volume = self.result[nout*self.q_input.nq + 1]/total_weight
        shell_volume = self.result[nout*self.q_input.nq + 2]/total_weight
        radius_effective = self.result[nout*self.q_input.nq + 3]/total_weight
        if shell_volume == 0.:
            shell_volume = 1.
        F1 = (self.result[1:nout*self.q_input.nq:nout]/total_weight
              if nout == 2 else None)
        F2 = self.result[0:nout*self.q_input.nq:nout]/total_weight
        return F1, F2, radius_effective, shell_volume, form_volume/shell_volume

    def release(self):
        # type: () -> None
        """
        Free resources associated with the kernel instance.
        """
        #print("null release kernel")
        pass

    def _call_kernel(self, call_details, values, cutoff, magnetic,
                     radius_effective_mode):
        # type: (CallDetails, np.ndarray, np.ndarray, float, bool, int) -> None
        """
        Call the kernel.  Subclasses defining kernels for particular execution
        engines need to provide an implementation for this.
        """
        raise NotImplementedError()
