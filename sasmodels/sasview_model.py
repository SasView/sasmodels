"""
Sasview model constructor.

Given a module defining an OpenCL kernel such as sasmodels.models.cylinder,
create a sasview model class to run that kernel as follows::

    from sasmodels.sasview_model import make_class
    from sasmodels.models import cylinder
    CylinderModel = make_class(cylinder, dtype='single')

The model parameters for sasmodels are different from those in sasview.
When reloading previously saved models, the parameters should be converted
using :func:`sasmodels.convert.convert`.
"""

# TODO: add a sasview=>sasmodels parameter translation layer
# this will allow us to use the new sasmodels as drop in replacements, and
# delay renaming parameters until all models have been converted.

import math
from copy import deepcopy
import warnings

import numpy as np

try:
    from .gpu import load_model
except ImportError,exc:
    warnings.warn(str(exc))
    warnings.warn("OpenCL not available --- using ctypes instead")
    from .dll import load_model


def make_class(kernel_module, dtype='single'):
    """
    Load the sasview model defined in *kernel_module*.

    Returns a class that can be used directly as a sasview model.
    """
    model =  load_model(kernel_module, dtype=dtype)
    def __init__(self, multfactor=1):
        SasviewModel.__init__(self, model)
    attrs = dict(__init__=__init__)
    ConstructedModel = type(model.info['name'],  (SasviewModel,), attrs)
    return ConstructedModel

class SasviewModel(object):
    """
    Sasview wrapper for opencl/ctypes model.
    """
    def __init__(self, model):
        """Initialization"""
        self._model = model

        self.name = model.info['name']
        self.description = model.info['description']
        self.category = None
        self.multiplicity_info = None
        self.is_multifunc = False

        ## interpret the parameters
        ## TODO: reorganize parameter handling
        self.details = dict()
        self.params = dict()
        self.dispersion = dict()
        partype = model.info['partype']
        for name,units,default,limits,ptype,description in model.info['parameters']:
            self.params[name] = default
            self.details[name] = [units]+limits

        for name in partype['pd-2d']:
            self.dispersion[name] = {
                'width': 0,
                'npts': 35,
                'nsigmas': 3,
                'type': 'gaussian',
            }

        self.orientation_params = (
            partype['orientation']
            + [n+'.width' for n in partype['orientation']]
            + partype['magnetic'])
        self.magnetic_params = partype['magnetic']
        self.fixed = [n+'.width' for n in partype['pd-2d']]
        self.non_fittable = []

        ## independent parameter name and unit [string]
        self.input_name = model.info.get("input_name","Q")
        self.input_unit = model.info.get("input_unit","A^{-1}")
        self.output_name = model.info.get("output_name","Intensity")
        self.output_unit = model.info.get("output_unit","cm^{-1}")

        ## _persistency_dict is used by sans.perspectives.fitting.basepage
        ## to store dispersity reference.
        ## TODO: _persistency_dict to persistency_dict throughout sasview
        self._persistency_dict = {}

        ## New fields introduced for opencl rewrite
        self.cutoff = 1e-5

    def __str__(self):
        """
        :return: string representation
        """
        return self.name

    def is_fittable(self, par_name):
        """
        Check if a given parameter is fittable or not

        :param par_name: the parameter name to check
        """
        return par_name.lower() in self.fixed
        #For the future
        #return self.params[str(par_name)].is_fittable()


    def getProfile(self):
        """
        Get SLD profile

        : return: (z, beta) where z is a list of depth of the transition points
                beta is a list of the corresponding SLD values
        """
        return None, None

    def setParam(self, name, value):
        """
        Set the value of a model parameter

        :param name: name of the parameter
        :param value: value of the parameter

        """
        # Look for dispersion parameters
        toks = name.split('.')
        if len(toks)==2:
            for item in self.dispersion.keys():
                if item.lower()==toks[0].lower():
                    for par in self.dispersion[item]:
                        if par.lower() == toks[1].lower():
                            self.dispersion[item][par] = value
                            return
        else:
            # Look for standard parameter
            for item in self.params.keys():
                if item.lower()==name.lower():
                    self.params[item] = value
                    return

        raise ValueError, "Model does not contain parameter %s" % name

    def getParam(self, name):
        """
        Set the value of a model parameter

        :param name: name of the parameter

        """
        # Look for dispersion parameters
        toks = name.split('.')
        if len(toks)==2:
            for item in self.dispersion.keys():
                if item.lower()==toks[0].lower():
                    for par in self.dispersion[item]:
                        if par.lower() == toks[1].lower():
                            return self.dispersion[item][par]
        else:
            # Look for standard parameter
            for item in self.params.keys():
                if item.lower()==name.lower():
                    return self.params[item]

        raise ValueError, "Model does not contain parameter %s" % name

    def getParamList(self):
        """
        Return a list of all available parameters for the model
        """
        list = self.params.keys()
        # WARNING: Extending the list with the dispersion parameters
        list.extend(self.getDispParamList())
        return list

    def getDispParamList(self):
        """
        Return a list of all available parameters for the model
        """
        # TODO: fix test so that parameter order doesn't matter
        ret = ['%s.%s'%(d.lower(), p)
               for d in self._model.info['partype']['pd-2d']
               for p in ('npts', 'nsigmas', 'width')]
        #print ret
        return ret

    def clone(self):
        """ Return a identical copy of self """
        return deepcopy(self)

    def run(self, x=0.0):
        """
        Evaluate the model

        :param x: input q, or [q,phi]

        :return: scattering function P(q)

        **DEPRECATED**: use calculate_Iq instead
        """
        if isinstance(x, (list,tuple)):
            q, phi = x
            return self.calculate_Iq([q * math.cos(phi)],
                                     [q * math.sin(phi)])[0]
        else:
            return self.calculate_Iq([float(x)])[0]


    def runXY(self, x=0.0):
        """
        Evaluate the model in cartesian coordinates

        :param x: input q, or [qx, qy]

        :return: scattering function P(q)

        **DEPRECATED**: use calculate_Iq instead
        """
        if isinstance(x, (list,tuple)):
            return self.calculate_Iq([float(x[0])],[float(x[1])])[0]
        else:
            return self.calculate_Iq([float(x)])[0]

    def evalDistribution(self, qdist):
        """
        Evaluate a distribution of q-values.

        * For 1D, a numpy array is expected as input: ::

            evalDistribution(q)

          where q is a numpy array.

        * For 2D, a list of numpy arrays are expected: [qx,qy],
          with 1D arrays::

              qx = [ qx[0], qx[1], qx[2], ....]

          and::

              qy = [ qy[0], qy[1], qy[2], ....]

        Then get ::

            q = numpy.sqrt(qx^2+qy^2)

        that is a qr in 1D array::

            q = [q[0], q[1], q[2], ....]


        :param qdist: ndarray of scalar q-values or list [qx,qy] where qx,qy are 1D ndarrays
        """
        if isinstance(qdist, (list,tuple)):
            # Check whether we have a list of ndarrays [qx,qy]
            qx, qy = qdist
            partype = self._model.info['partype']
            if not partype['orientation'] and not partype['magnetic']:
                return self.calculate_Iq(np.sqrt(qx**2+qy**2))
            else:
                return self.calculate_Iq(qx, qy)

        elif isinstance(qdist, np.ndarray):
            # We have a simple 1D distribution of q-values
            return self.calculate_Iq(qdist)

        else:
            raise TypeError("evalDistribution expects q or [qx, qy], not %r"%type(qdist))

    def calculate_Iq(self, *args):
        """
        Calculate Iq for one set of q with the current parameters.

        If the model is 1D, use *q*.  If 2D, use *qx*, *qy*.

        This should NOT be used for fitting since it copies the *q* vectors
        to the card for each evaluation.
        """
        q_vectors = [np.asarray(q) for q in args]
        fn = self._model(self._model.make_input(q_vectors))
        pars = [self.params[v] for v in fn.fixed_pars]
        pd_pars = [self._get_weights(p) for p in fn.pd_pars]
        result = fn(pars, pd_pars, self.cutoff)
        fn.input.release()
        fn.release()
        return result

    def calculate_ER(self):
        """
        Calculate the effective radius for P(q)*S(q)

        :return: the value of the effective radius
        """
        ER = self._model.info.get('ER', None)
        if ER is None:
            return 1.0
        else:
            vol_pars = self._model.info['partype']['volume']
            values, weights = self._dispersion_mesh(vol_pars)
            fv = ER(*values)
            #print values[0].shape, weights.shape, fv.shape
            return np.sum(weights*fv) / np.sum(weights)

    def calculate_VR(self):
        """
        Calculate the volf ratio for P(q)*S(q)

        :return: the value of the volf ratio
        """
        VR = self._model.info.get('VR', None)
        if VR is None:
            return 1.0
        else:
            vol_pars = self._model.info['partype']['volume']
            values, weights = self._dispersion_mesh(vol_pars)
            whole,part = VR(*values)
            return np.sum(weights*part)/np.sum(weights*whole)

    def set_dispersion(self, parameter, dispersion):
        """
        Set the dispersion object for a model parameter

        :param parameter: name of the parameter [string]
        :param dispersion: dispersion object of type Dispersion
        """
        if parameter.lower() in (s.lower() for s in self.params.keys()):
            # TODO: Store the disperser object directly in the model.
            # The current method of creating one on the fly whenever it is
            # needed is kind of funky.
            # Note: can't seem to get disperser parameters from sasview
            # (1) Could create a sasview model that has not yet # been
            # converted, assign the disperser to one of its polydisperse
            # parameters, then retrieve the disperser parameters from the
            # sasview model.  (2) Could write a disperser parameter retriever
            # in sasview.  (3) Could modify sasview to use sasmodels.weights
            # dispersers.
            # For now, rely on the fact that the sasview only ever uses
            # new dispersers in the set_dispersion call and create a new
            # one instead of trying to assign parameters.
            from . import weights
            disperser = weights.dispersers[dispersion.__class__.__name__]
            dispersion = weights.models[disperser]()
            self.dispersion[parameter] = dispersion.get_pars()
        else:
            raise ValueError("%r is not a dispersity or orientation parameter")

    def _dispersion_mesh(self, pars):
        """
        Create a mesh grid of dispersion parameters and weights.

        Returns [p1,p2,...],w where pj is a vector of values for parameter j
        and w is a vector containing the products for weights for each
        parameter set in the vector.
        """
        values, weights = zip(*[self._get_weights(p) for p in pars])
        values = [v.flatten() for v in np.meshgrid(*values)]
        weights = np.vstack([v.flatten() for v in np.meshgrid(*weights)])
        weights = np.prod(weights, axis=0)
        return values, weights

    def _get_weights(self, par):
        from . import weights

        relative = self._model.info['partype']['pd-rel']
        limits = self._model.info['limits']
        dis = self.dispersion[par]
        v,w = weights.get_weights(
            dis['type'], dis['npts'], dis['width'], dis['nsigmas'],
            self.params[par], limits[par], par in relative)
        return v,w/w.max()

