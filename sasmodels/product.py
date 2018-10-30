"""
Product model
-------------

The product model multiplies the structure factor by the form factor,
modulated by the effective radius of the form.  The resulting model
has a attributes of both the model description (with parameters, etc.)
and the module evaluator (with call, release, etc.).

To use it, first load form factor P and structure factor S, then create
*make_product_info(P, S)*.
"""
from __future__ import print_function, division

from collections import OrderedDict

from copy import copy
import numpy as np  # type: ignore

from .modelinfo import ParameterTable, ModelInfo, Parameter, parse_parameter
from .kernel import KernelModel, Kernel
from .details import make_details, dispersion_mesh

# pylint: disable=unused-import
try:
    from typing import Tuple, Callable, Union
except ImportError:
    pass
else:
    from .modelinfo import ParameterSet
# pylint: enable=unused-import

# TODO: make estimates available to constraints
#ESTIMATED_PARAMETERS = [
#    ["est_radius_effective", "A", 0.0, [0, np.inf], "", "Estimated effective radius"],
#    ["est_volume_ratio", "", 1.0, [0, np.inf], "", "Estimated volume ratio"],
#]

STRUCTURE_MODE_ID = "structure_factor_mode"
RADIUS_MODE_ID = "radius_effective_mode"
RADIUS_ID = "radius_effective"
VOLFRAC_ID = "volfraction"
def make_extra_pars(p_info):
    pars = []
    if p_info.have_Fq:
        par = parse_parameter(
                STRUCTURE_MODE_ID,
                "",
                0,
                [["P*S","P*(1+beta*(S-1))"]],
                "",
                "Structure factor calculation")
        pars.append(par)
    if p_info.effective_radius_type is not None:
        par = parse_parameter(
                RADIUS_MODE_ID,
                "",
                0,
                [["unconstrained"] + p_info.effective_radius_type],
                "",
                "Effective radius calculation")
        pars.append(par)
    return pars

def make_product_info(p_info, s_info):
    # type: (ModelInfo, ModelInfo) -> ModelInfo
    """
    Create info block for product model.
    """
    # Make sure effective radius is the first parameter and
    # make sure volume fraction is the second parameter of the
    # structure factor calculator.  Structure factors should not
    # have any magnetic parameters
    if not len(s_info.parameters.kernel_parameters) >= 2:
        raise TypeError("S needs {} and {} as its first parameters".format(RADIUS_ID, VOLFRAC_ID))
    if not s_info.parameters.kernel_parameters[0].id == RADIUS_ID:
        raise TypeError("S needs {} as first parameter".format(RADIUS_ID))
    if not s_info.parameters.kernel_parameters[1].id == VOLFRAC_ID:
        raise TypeError("S needs {} as second parameter".format(VOLFRAC_ID))
    if not s_info.parameters.magnetism_index == []:
        raise TypeError("S should not have SLD parameters")
    p_id, p_name, p_pars = p_info.id, p_info.name, p_info.parameters
    s_id, s_name, s_pars = s_info.id, s_info.name, s_info.parameters

    # Create list of parameters for the combined model.  If there
    # are any names in P that overlap with those in S, modify the name in S
    # to distinguish it.
    p_set = set(p.id for p in p_pars.kernel_parameters)
    s_list = [(_tag_parameter(par) if par.id in p_set else par)
              for par in s_pars.kernel_parameters]
    # Check if still a collision after renaming.  This could happen if for
    # example S has volfrac and P has both volfrac and volfrac_S.
    if any(p.id in p_set for p in s_list):
        raise TypeError("name collision: P has P.name and P.name_S while S has S.name")

    # make sure effective radius is not a polydisperse parameter in product
    s_list[0] = copy(s_list[0])
    s_list[0].polydisperse = False

    translate_name = dict((old.id, new.id) for old, new
                          in zip(s_pars.kernel_parameters, s_list))
    combined_pars = p_pars.kernel_parameters + s_list + make_extra_pars(p_info)
    parameters = ParameterTable(combined_pars)
    parameters.max_pd = p_pars.max_pd + s_pars.max_pd
    def random():
        combined_pars = p_info.random()
        s_names = set(par.id for par in s_pars.kernel_parameters)
        combined_pars.update((translate_name[k], v)
                             for k, v in s_info.random().items()
                             if k in s_names)
        return combined_pars

    model_info = ModelInfo()
    model_info.id = '@'.join((p_id, s_id))
    model_info.name = '@'.join((p_name, s_name))
    model_info.filename = None
    model_info.title = 'Product of %s and %s'%(p_name, s_name)
    model_info.description = model_info.title
    model_info.docs = model_info.title
    model_info.category = "custom"
    model_info.parameters = parameters
    model_info.random = random
    #model_info.single = p_info.single and s_info.single
    model_info.structure_factor = False
    model_info.variant_info = None
    #model_info.tests = []
    #model_info.source = []
    # Remember the component info blocks so we can build the model
    model_info.composition = ('product', [p_info, s_info])
    model_info.control = p_info.control
    model_info.hidden = p_info.hidden
    if getattr(p_info, 'profile', None) is not None:
        profile_pars = set(p.id for p in p_info.parameters.kernel_parameters)
        def profile(**kwargs):
            # extract the profile args
            kwargs = dict((k, v) for k, v in kwargs.items() if k in profile_pars)
            return p_info.profile(**kwargs)
    else:
        profile = None
    model_info.profile = profile
    model_info.profile_axes = p_info.profile_axes

    # TODO: delegate random to p_info, s_info
    #model_info.random = lambda: {}

    ## Show the parameter table
    #from .compare import get_pars, parlist
    #print("==== %s ====="%model_info.name)
    #values = get_pars(model_info)
    #print(parlist(model_info, values, is2d=True))
    return model_info

def _tag_parameter(par):
    """
    Tag the parameter name with _S to indicate that the parameter comes from
    the structure factor parameter set.  This is only necessary if the
    form factor model includes a parameter of the same name as a parameter
    in the structure factor.
    """
    par = copy(par)
    # Protect against a vector parameter in S by appending the vector length
    # to the renamed parameter.  Note: haven't tested this since no existing
    # structure factor models contain vector parameters.
    vector_length = par.name[len(par.id):]
    par.id = par.id + "_S"
    par.name = par.id + vector_length
    return par

def _intermediates(
        F1,               # type: np.ndarray
        F2,               # type: np.ndarray
        S,                # type: np.ndarray
        scale,            # type: float
        effective_radius, # type: float
        beta_mode,        # type: bool
        ):
    # type: (...) -> OrderedDict[str, Union[np.ndarray, float]]
    """
    Returns intermediate results for beta approximation-enabled product.
    The result may be an array or a float.
    """
    if beta_mode:
        # TODO: 1. include calculated Q vector
        # TODO: 2. consider implications if there are intermediate results in P(Q)
        parts = OrderedDict((
            ("P(Q)", scale*F2),
            ("S(Q)", S),
            ("beta(Q)", F1**2 / F2),
            ("S_eff(Q)", 1 + (F1**2 / F2)*(S-1)),
            ("effective_radius", effective_radius),
            # ("I(Q)", scale*(F2 + (F1**2)*(S-1)) + bg),
        ))
    else:
        parts = OrderedDict((
            ("P(Q)", scale*F2),
            ("S(Q)", S),
            ("effective_radius", effective_radius),
        ))
    return parts

class ProductModel(KernelModel):
    def __init__(self, model_info, P, S):
        # type: (ModelInfo, KernelModel, KernelModel) -> None
        #: Combined info plock for the product model
        self.info = model_info
        #: Form factor modelling individual particles.
        self.P = P
        #: Structure factor modelling interaction between particles.
        self.S = S

        #: Model precision. This is not really relevant, since it is the
        #: individual P and S models that control the effective dtype,
        #: converting the q-vectors to the correct type when the kernels
        #: for each are created. Ideally this should be set to the more
        #: precise type to avoid loss of precision, but precision in q is
        #: not critical (single is good enough for our purposes), so it just
        #: uses the precision of the form factor.
        self.dtype = P.dtype  # type: np.dtype

    def make_kernel(self, q_vectors):
        # type: (List[np.ndarray]) -> Kernel
        # Note: may be sending the q_vectors to the GPU twice even though they
        # are only needed once.  It would mess up modularity quite a bit to
        # handle this optimally, especially since there are many cases where
        # separate q vectors are needed (e.g., form in python and structure
        # in opencl; or both in opencl, but one in single precision and the
        # other in double precision).

        p_kernel = self.P.make_kernel(q_vectors)
        s_kernel = self.S.make_kernel(q_vectors)
        return ProductKernel(self.info, p_kernel, s_kernel)

    def release(self):
        # type: (None) -> None
        """
        Free resources associated with the model.
        """
        self.P.release()
        self.S.release()


class ProductKernel(Kernel):
    def __init__(self, model_info, p_kernel, s_kernel):
        # type: (ModelInfo, Kernel, Kernel) -> None
        self.info = model_info
        self.p_kernel = p_kernel
        self.s_kernel = s_kernel
        self.dtype = p_kernel.dtype
        self.results = []  # type: List[np.ndarray]

    def __call__(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, float, bool) -> np.ndarray

        p_info, s_info = self.info.composition[1]
        p_npars = p_info.parameters.npars
        p_length = call_details.length[:p_npars]
        p_offset = call_details.offset[:p_npars]
        s_npars = s_info.parameters.npars
        s_length = call_details.length[p_npars:p_npars+s_npars]
        s_offset = call_details.offset[p_npars:p_npars+s_npars]

        # Beta mode parameter is the first parameter after P and S parameters
        have_beta_mode = p_info.have_Fq
        beta_mode_offset = 2+p_npars+s_npars
        beta_mode = (values[beta_mode_offset] > 0) if have_beta_mode else False
        if beta_mode and self.p_kernel.dim== '2d':
            raise NotImplementedError("beta not yet supported for 2D")

        # R_eff type parameter is the second parameter after P and S parameters
        # unless the model doesn't support beta mode, in which case it is first
        have_radius_type = p_info.effective_radius_type is not None
        radius_type_offset = 2+p_npars+s_npars + (1 if have_beta_mode else 0)
        radius_type = int(values[radius_type_offset]) if have_radius_type else 0

        # Retrieve the volume fraction, which is the second of the
        # 'S' parameters in the parameter list, or 2+np in 0-origin,
        # as well as the scale and background.
        volfrac = values[3+p_npars]
        scale, background = values[0], values[1]

        # if there are magnetic parameters, they will only be on the
        # form factor P, not the structure factor S.
        nmagnetic = len(self.info.parameters.magnetism_index)
        if nmagnetic:
            spin_index = self.info.parameters.npars + 2
            magnetism = values[spin_index: spin_index+3+3*nmagnetic]
        else:
            magnetism = []
        nvalues = self.info.parameters.nvalues
        nweights = call_details.num_weights
        weights = values[nvalues:nvalues + 2*nweights]

        # Construct the calling parameters for P.
        p_details = make_details(p_info, p_length, p_offset, nweights)
        p_values = [
            [1., 0.], # scale=1, background=0,
            values[2:2+p_npars],
            magnetism,
            weights]
        spacer = (32 - sum(len(v) for v in p_values)%32)%32
        p_values.append([0.]*spacer)
        p_values = np.hstack(p_values).astype(self.p_kernel.dtype)

        # Construct the calling parameters for S.
        if radius_type > 0:
            # If R_eff comes from form factor, make sure it is monodisperse.
            # weight is set to 1 later, after the value array is created
            s_length[0] = 1
        s_details = make_details(s_info, s_length, s_offset, nweights)
        s_values = [
            [1., 0.], # scale=1, background=0,
            values[2+p_npars:2+p_npars+s_npars],
            weights,
        ]
        spacer = (32 - sum(len(v) for v in s_values)%32)%32
        s_values.append([0.]*spacer)
        s_values = np.hstack(s_values).astype(self.s_kernel.dtype)

        # Call the form factor kernel to compute <F> and <F^2>.
        # If the model doesn't support Fq the returned <F> will be None.
        F1, F2, effective_radius, shell_volume, volume_ratio = self.p_kernel.Fq(
            p_details, p_values, cutoff, magnetic, radius_type)

        # Call the structure factor kernel to compute S.
        # Plug R_eff from the form factor into structure factor parameters
        # and scale volume fraction by form:shell volume ratio. These changes
        # needs to be both in the initial value slot as well as the
        # polydispersity distribution slot in the values array due to
        # implementation details in kernel_iq.c.
        #print("R_eff=%d:%g, volfrac=%g, volume ratio=%g"%(radius_type, effective_radius, volfrac, volume_ratio))
        if radius_type > 0:
            # set the value to the model R_eff and set the weight to 1
            s_values[2] = s_values[2+s_npars+s_offset[0]] = effective_radius
            s_values[2+s_npars+s_offset[0]+nweights] = 1.0
        s_values[3] = s_values[2+s_npars+s_offset[1]] = volfrac*volume_ratio
        S = self.s_kernel.Iq(s_details, s_values, cutoff, False)

        # Determine overall scale factor. Hollow shapes are weighted by
        # shell_volume, so that is needed for volume normalization.  For
        # solid shapes we can use shell_volume as well since it is equal
        # to form volume.
        combined_scale = scale*volfrac/shell_volume

        # Combine form factor and structure factor
        #print("beta", beta_mode, F1, F2, S)
        PS = F2 + F1**2*(S-1) if beta_mode else F2*S
        final_result = combined_scale*PS + background

        # Capture intermediate values so user can see them.  These are
        # returned as a lazy evaluator since they are only needed in the
        # GUI, and not for each evaluation during a fit.
        # TODO: return the results structure with the final results
        # That way the model calcs are idempotent. Further, we can
        # generalize intermediates to various other model types if we put it
        # kernel calling interface.  Could do this as an "optional"
        # return value in the caller, though in that case we could return
        # the results directly rather than through a lazy evaluator.
        self.results = lambda: _intermediates(
            F1, F2, S, combined_scale, effective_radius, beta_mode)

        return final_result

    def release(self):
        # type: () -> None
        self.p_kernel.release()
        self.s_kernel.release()
