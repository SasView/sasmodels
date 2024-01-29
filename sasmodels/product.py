"""
Product model
-------------

The product model multiplies the structure factor by the form factor,
modulated by the effective radius of the form.  The resulting model
has a attributes of both the model description (with parameters, etc.)
and the module evaluator (with call, release, etc.).

To use it, first load form factor P and structure factor S, then create
*make_product_info(P, S)*.

The P@S models is somewhat complicated because there are many special
parameters that need to be handled in particular ways.  Much of the
code is used to figure out what special parameters we have, where to
find them in the P@S model inputs and how to distribute them to the underlying
P and S model calculators.

The parameter packet received by the P@S is a :class:`.details.CallDetails`
structure along with a data vector. The CallDetails structure indicates which
parameters are polydisperse, the length of the distribution, and where to
find it in the data vector.  The distributions are ordered from longest to
shortest, with length 1 distributions filling out the distribution set.  That
way the kernel calculator doesn't have to check if it needs another nesting
level since it is always there.  The data vector consists of a list of target
values for the parameters, followed by a concatenation of the distribution
values, and then followed by a concatenation of the distribution weights.
Given the combined details and data for P@S, we must decompose them in to
details for P and details for S separately, which unfortunately requires
intimate knowledge of the data structures and tricky code.

The special parameters are:

* *scale* and *background*:
    First two parameters of the value list in each of P, S and P@S.
    When decomposing P@S parameters, ignore *scale* and *background*,
    instead using 1 and 0 for the first two slots of both P and S.
    After calling P and S individually, the results are combined as
    :code:`volfraction*scale*P*S + background`.  The *scale* and
    *background* do not show up in the polydispersity structure so
    they are easy to handle.

* *volfraction*:
    Always the first parameter of S, but it may also be in P. If it is in P,
    then *P.volfraction* is used in the combined P@S list, and
    *S.volfraction* is elided, otherwise *S.volfraction* is used. If we
    are using *volfraction* from P we can treat it like all the other P
    parameters when calling P, but when calling S we need to insert the
    *P.volfraction* into data vector for S and assign a slot of length 1
    in the distribution. Because we are using the original layout of the
    distribution vectors from P@S, but copying it into private data
    vectors for S and P, we are free to "borrow" a P slots to store the
    missing *S.volfraction* distribution.  We use the *P.volfraction*
    slot itself but any slot will work.

    For hollow shapes, *volfraction* represents the volume fraction of
    material but S needs the volume fraction enclosed by the shape. The
    answer is to scale the user specified volume fraction by the form:shell
    ratio computed from the average form volume and average shell volume
    returned from P. Use the original *volfraction* divided by *shell_volume*
    to compute the number density, and scale P@S by that to get absolute
    scaling on the final *I(q)*. The *scale* for P@S should therefore usually
    be one.

* *radius_effective*:
    Always the second parameter of S and always part of P@S, but never in P.
    The value may be calculated using *P.radius_effective()* or it
    may be set to the *radius_effective* value in P@S, depending on
    *radius_effective_mode*.  If part of S, the value may be polydisperse.
    If calculated by P, then it will be the weighted average of effective
    radii computed for the polydisperse shape parameters.

* *structure_factor_mode*
    If P@S supports beta approximation (i.e., if it has the *Fq* function that
    returns <FF*> and <F><F*>), then *structure_factor_mode* will be added
    to the P@S parameters right after the S parameters.  This mode may be 0
    for the monodisperse approximation or 1 for the beta approximation.  We
    will add more values here as we implemented more complicated operations,
    but for now P and S must be computed separately.  If beta, then we
    return *I = scale volfrac/volume ( <FF> + <F>^2 (S-1)) + background*.
    If not beta then return *I = scale/volume P S + background* .  In both
    cases, return the appropriate immediate values.

* *radius_effective_mode*
    If P defines the *radius_effective* function (and therefore
    *P.info.radius_effective_modes* is a list of effective radius modes),
    then *radius_effective_mode* will be the final parameter in P@S.  Mode
    will be zero if *radius_effective* is defined by the user using the S
    parameter; any other value and the *radius_effective* parameter will be
    filled in from the value computed in P.  In the latter case, the
    polydispersity information for *S.radius_effective* will need to be
    suppressed, with pd length set to 1, the first value set to the
    effective radius and the first weight set to 1.  Do this after composing
    the S data vector so the inputs are left untouched.

* *regular parameters*
    The regular P parameters form a block of length *P.info.npars* at the
    start of the data vector (after scale and background).  These will be
    followed by *S.effective_radius*, and *S.volfraction* (if *P.volfraction*
    is absent), and then the regular S parameters.  The P and S blocks can
    be copied as a group into the respective P and S data vectors.
    We can copy the distribution value and weight vectors untouched to both
    the P and S data vectors since they are referenced by offset and length.
    We can update the radius_effective slots in the P data vector with
    *P.radius_effective()* if needed.

* *magnetic parameters*
    For each P parameter that is an SLD there will be a set of three magnetic
    parameters tacked on to P@S after the regular P and S and after the
    special *structure_factor_mode* and *radius_effective_mode*.  These
    can be copied as a group after the regular P parameters.  There won't
    be any magnetic S parameters.

"""
from __future__ import print_function, division

from collections import OrderedDict

from copy import copy
import numpy as np  # type: ignore

from .modelinfo import ParameterTable, ModelInfo, parse_parameter
from .modelinfo import NUM_MAGFIELD_PARS, NUM_MAGNETIC_PARS, NUM_COMMON_PARS
from .kernel import KernelModel, Kernel
from .details import make_details

# pylint: disable=unused-import
try:
    from typing import OrderedDict as OrderedDictType
    import typing
    from typing import Tuple, Callable, Union, List, Optional, Dict
    from .modelinfo import ParameterSet, Parameter
    from .details import CallDetails
    Parts = Dict[str, Union[float, np.ndarray, Tuple[np.ndarray, np.ndarray]]]
except ImportError:
    pass
# pylint: enable=unused-import

# TODO: make shape averages available to constraints
#ESTIMATED_PARAMETERS = [
#    ["mean_radius_effective", "A", 0.0, [0, np.inf], "", "mean effective radius"],
#    ["mean_volume", "A", 0.0, [0, np.inf], "", "mean volume"],
#    ["mean_volume_ratio", "", 1.0, [0, np.inf], "", "mean form: mean shell volume ratio"],
#]
STRUCTURE_MODE_ID = "structure_factor_mode"
RADIUS_MODE_ID = "radius_effective_mode"
RADIUS_ID = "radius_effective"
VOLFRAC_ID = "volfraction"
def make_extra_pars(p_info):
    # type: (ModelInfo) -> List[Parameter]
    """
    Create parameters for structure factor and effective radius modes.
    """
    pars = []
    if p_info.have_Fq:
        par = parse_parameter(
            STRUCTURE_MODE_ID,
            "",
            0,
            [["P*S", "P*(1+beta*(S-1))"]],
            "",
            "Structure factor calculation")
        pars.append(par)
    if p_info.radius_effective_modes is not None:
        par = parse_parameter(
            RADIUS_MODE_ID,
            "",
            1,
            [["unconstrained"] + p_info.radius_effective_modes],
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
    if RADIUS_ID in p_info.parameters:
        raise TypeError("P should not have {}".format(RADIUS_ID))
    p_id, p_name, p_pars = p_info.id, p_info.name, p_info.parameters
    s_id, s_name, s_pars = s_info.id, s_info.name, s_info.parameters
    p_has_volfrac = VOLFRAC_ID in p_info.parameters

    # Create list of parameters for the combined model.  If a name in
    # S overlaps a name in P, tag the S parameter name to distinguish it.
    # If the tagged name also collides it will be caught by the parameter
    # table builder.  Similarly if any special names are abused.  Need the
    # pairs to create the translation table for random model generation.
    p_set = set(p.id for p in p_pars.kernel_parameters)
    s_pairs = [(par, (_tag_parameter(par) if par.id in p_set else par))
               for par in s_pars.kernel_parameters
               # Note: exclude volfraction from s_list if volfraction in p
               if par.id != VOLFRAC_ID or not p_has_volfrac]
    s_list = [pair[0] for pair in s_pairs]

    # Build combined parameter table
    combined_pars = p_pars.kernel_parameters + s_list + make_extra_pars(p_info)
    parameters = ParameterTable(combined_pars)
    # Allow for the scenario in which each component has all its PD parameters
    # active simultaneously.  details.make_details() will throw an error if
    # too many are used from any one component.
    parameters.max_pd = p_pars.max_pd + s_pars.max_pd

    # TODO: does user-defined polydisperse S.radius_effective make sense?
    # make sure effective radius is not a polydisperse parameter in product
    #s_list[0] = copy(s_list[0])
    #s_list[0].polydisperse = False

    s_translate = {old.id: new.id for old, new in s_pairs}
    def random():
        """Random set of model parameters for product model"""
        combined_pars = p_info.random()
        combined_pars.update((s_translate[k], v)
                             for k, v in s_info.random().items()
                             if k in s_translate)
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
    #model_info.tests = []
    #model_info.source = []
    # Remember the component info blocks so we can build the model
    model_info.composition = ('product', [p_info, s_info])
    model_info.hidden = p_info.hidden
    if getattr(p_info, 'profile', None) is not None:
        profile_pars = set(p.id for p in p_info.parameters.kernel_parameters)
        def profile(**kwargs):
            """Return SLD profile of the form factor as a function of radius."""
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

def _intermediates(Q, F, Fsq, S, scale, volume, volume_ratio, radius_effective,
                   beta_mode, P_intermediate):
    # type: (np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, float, float, bool, Optiona[Callable[[], Parts]]) -> Parts
    """
    Returns intermediate results for beta approximation-enabled product.
    The result may be an array or a float.
    """
    parts = OrderedDict()  # type: Parts
    parts["P(Q)"] = (Q, scale*Fsq)
    if P_intermediate is not None:
        parts["P(Q) parts"] = P_intermediate()
    parts["volume"] = volume
    parts["volume_ratio"] = volume_ratio
    parts["radius_effective"] = radius_effective
    parts["S(Q)"] = (Q, S)
    if beta_mode:
        parts["beta(Q)"] = (Q, F**2 / Fsq)
        parts["S_eff(Q)"] = (Q, 1 + (F**2 / Fsq)*(S-1))
        #parts["I(Q)", scale*(Fsq + (F**2)*(S-1)) + bg
    return parts

class ProductModel(KernelModel):
    """
    Model definition for product model.
    """
    info = None  # type: ModelInfo
    P = None  # type: KernelModel
    S = None  # type: KernelModel
    dtype = None  # type: np.dtype
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
        return ProductKernel(self.info, p_kernel, s_kernel, q_vectors)
    make_kernel.__doc__ = KernelModel.make_kernel.__doc__

    def release(self):
        # type: (None) -> None
        """
        Free resources associated with the model.
        """
        self.P.release()
        self.S.release()


class ProductKernel(Kernel):
    """
    Instantiated kernel for product model.
    """
    def __init__(self, model_info, p_kernel, s_kernel, q):
        # type: (ModelInfo, Kernel, Kernel, Tuple[np.ndarray]) -> None
        self.info = model_info
        self.q = q
        self.p_kernel = p_kernel
        self.s_kernel = s_kernel
        self.dtype = p_kernel.dtype
        self.results = None  # type: Callable[[], OrderedDict]

        # Find index of volfraction parameter in parameter list
        for k, p in enumerate(model_info.parameters.call_parameters):
            if p.id == VOLFRAC_ID:
                self._volfrac_index = k
                break
        else:
            raise RuntimeError("no %s parameter in %s"%(VOLFRAC_ID, self))

        p_info, s_info = self.info.composition[1]
        p_npars = p_info.parameters.npars
        s_npars = s_info.parameters.npars

        have_beta_mode = p_info.have_Fq
        have_er_mode = p_info.radius_effective_modes is not None
        volfrac_in_p = self._volfrac_index < p_npars + NUM_COMMON_PARS

        # Slices into the details length/offset structure for P@S.
        # Made complicated by the possibly missing volfraction in S.
        self._p_detail_slice = slice(0, p_npars)
        self._s_detail_slice = slice(p_npars, p_npars+s_npars-volfrac_in_p)
        self._volfrac_in_p = volfrac_in_p

        # P block from data vector, without scale and background
        first_p = NUM_COMMON_PARS
        last_p = p_npars + NUM_COMMON_PARS
        self._p_value_slice = slice(first_p, last_p)

        # radius_effective is the first parameter in S from the data vector.
        self._er_index = last_p

        # S block from data vector, without scale, background, volfrac or er.
        first_s = last_p + NUM_COMMON_PARS - volfrac_in_p
        last_s = first_s + s_npars - NUM_COMMON_PARS
        self._s_value_slice = slice(first_s, last_s)

        # S distribution block in S data vector starts after all S values
        self._s_dist_slice = slice(NUM_COMMON_PARS + s_npars, None)

        # structure_factor_mode is the first parameter after P and S.  Skip
        # 2 for scale and background, and subtract 1 in case there is no
        # volfraction in S.
        self._beta_mode_index = last_s if have_beta_mode else 0

        # radius_effective_mode is the second parameter after P and S
        # unless structure_factor_mode isn't available, in which case it
        # is first.
        self._er_mode_index = last_s + have_beta_mode if have_er_mode else 0

        # Magnetic parameters are after everything else.  If they exist,
        # they will only be for form factor P, not structure factor S.
        first_mag = last_s + have_beta_mode + have_er_mode
        mag_pars = NUM_MAGNETIC_PARS*p_info.parameters.nmagnetic
        last_mag = first_mag + (mag_pars + NUM_MAGFIELD_PARS if mag_pars else 0)
        self._magentic_slice = slice(first_mag, last_mag)

    def Iq(self, call_details, values, cutoff, magnetic):
        # type: (CallDetails, np.ndarray, float, bool) -> np.ndarray
        p_info, s_info = self.info.composition[1]

        # Retrieve values from the data vector
        scale, background = values[0], values[1]
        volfrac = values[self._volfrac_index]
        er_mode = (int(values[self._er_mode_index])
                   if self._er_mode_index > 0 else 0)
        beta_mode = (values[self._beta_mode_index] > 0
                     if self._beta_mode_index > 0 else False)

        nvalues = self.info.parameters.nvalues
        nweights = call_details.num_weights
        weights = values[nvalues:nvalues + 2*nweights]

        # Can't do 2d and beta_mode just yet
        if beta_mode and self.p_kernel.dim == '2d':
            raise NotImplementedError("beta not yet supported for 2D")

        # Construct the calling parameters for P.
        p_length = call_details.length[self._p_detail_slice]
        p_offset = call_details.offset[self._p_detail_slice]
        p_details = make_details(p_info, p_length, p_offset, nweights)
        p_values = [
            [1., 0.], # scale=1, background=0,
            values[self._p_value_slice],
            values[self._magentic_slice],
            weights]
        spacer = (32 - sum(len(v) for v in p_values)%32)%32
        p_values.append([0.]*spacer)
        p_values = np.hstack(p_values).astype(self.p_kernel.dtype)

        # Call the form factor kernel to compute <F> and <F^2>.
        # If the model doesn't support Fq the returned <F> will be None.
        F, Fsq, radius_effective, shell_volume, volume_ratio \
            = self.p_kernel.Fq(p_details, p_values, cutoff, magnetic, er_mode)
        p_intermediate = getattr(self.p_kernel, 'results', None)

        # TODO: async call to the GPU

        # Construct the calling parameters for S.
        s_length = call_details.length[self._s_detail_slice]
        s_offset = call_details.offset[self._s_detail_slice]
        if self._volfrac_in_p:
            # Volfrac is in P and missing from S so insert a slot for it.  Say
            # the distribution is length 1 and use the slot for volfraction
            # from the P distribution.
            s_length = np.insert(s_length, 1, 1)
            s_offset = np.insert(s_offset, 1, p_offset[self._volfrac_index - 2])
        if er_mode > 0:
            # If effective_radius comes from P, make sure it is monodisperse.
            # Weight is set to 1 later, after the value array is created
            s_length[0] = 1
        s_details = make_details(s_info, s_length, s_offset, nweights)
        s_values = [
            [1., # scale=1
             0., # background=0,
             values[self._er_index], # S.radius_effective; may be replaced by P
             0.], # volfraction; will be replaced by volfrac * volume_ratio
            # followed by S parameters after effective_radius and volfraction
            values[self._s_value_slice],
            weights,
        ]
        spacer = (32 - sum(len(v) for v in s_values)%32)%32
        s_values.append([0.]*spacer)
        s_values = np.hstack(s_values).astype(self.s_kernel.dtype)

        # Plug R_eff from the form factor into structure factor parameters
        # and scale volume fraction by form:shell volume ratio. These changes
        # needs to be both in the initial value slot as well as the
        # polydispersity distribution slot in the values array due to
        # implementation details in kernel_iq.c.
        #print("R_eff=%d:%g, volfrac=%g, volume ratio=%g"
        #      % (radius_type, radius_effective, volfrac, volume_ratio))
        s_dist = s_values[self._s_dist_slice]
        if er_mode > 0:
            # set the value to the model R_eff and set the weight to 1
            s_values[NUM_COMMON_PARS] = s_dist[s_offset[0]] = radius_effective
            s_dist[s_offset[0]+nweights] = 1.0
        s_values[NUM_COMMON_PARS+1] = s_dist[s_offset[1]] = volfrac*volume_ratio
        s_dist[s_offset[1]+nweights] = 1.0

        # Call the structure factor kernel to compute S.
        S = self.s_kernel.Iq(s_details, s_values, cutoff, False)
        #print("P", Fsq[:10])
        #print("S", S[:10])
        #print(radius_effective, volfrac*volume_ratio)

        # Combine form factor and structure factor
        #print("beta", beta_mode, F, Fsq, S)
        PS = Fsq + F**2*(S-1) if beta_mode else Fsq*S

        # Determine overall scale factor. Hollow shapes are weighted by
        # shell_volume, so that is needed for number density estimation.
        # For solid shapes we can use shell_volume as well since it is
        # equal to form volume.  If P already has a volfraction parameter,
        # then assume that it is already on absolute scale, and don't
        # include volfrac in the combined_scale.
        combined_scale = scale/shell_volume
        if not self._volfrac_in_p:
            combined_scale *= volfrac
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
            self.q, F, Fsq, S, combined_scale, shell_volume, volume_ratio,
            radius_effective, beta_mode, p_intermediate)

        return final_result

    Iq.__doc__ = Kernel.Iq.__doc__
    __call__ = Iq

    def release(self):
        # type: () -> None
        """Free resources associated with the kernel."""
        self.p_kernel.release()
        self.s_kernel.release()
