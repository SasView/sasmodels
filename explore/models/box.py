"""
Pure python parallelepiped (2d only)
"""
from numpy import inf

from sasmodels.special import sas_sinx_x

name = "box"
title = "pure python box"
description = """
        Model for prisms of different cross-sections with orientation average using the Fibonacci quadrature"""
category = "shape:cylinder"
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Parallelepiped scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["length_a", "Ang", 35, [0, inf], "volume",
               "Shorter side of the parallelepiped"],
              ["length_b", "Ang", 75, [0, inf], "volume",
               "Second side of the parallelepiped"],
              ["length_c", "Ang", 400, [0, inf], "volume",
               "Larger side of the parallelepiped"],
              ["theta", "degrees", 60, [-180, 180], "orientation",
               "c axis to beam angle"],
              ["phi", "degrees", 60, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 60, [-360, 360], "orientation",
               "rotation about c axis"],
             ]

### Functions for geometrical calculations: volume, surface of the cross section, average radius, edge length etc.
def form_volume(length_a, length_b, length_c):
    """
    Computes the volume of a parallelepiped.

    Returns
    -------
    volume : float
    """
    return length_a*length_b*length_c

def _fq(qa, qb, qc, length_a, length_b, length_c):
    """
    Computes the form factor amplitude of a cylinder

    Parameters
    ----------
    qa, qb, qc : np.ndarray shape (N,)
        components of the scattering vector q

    Returns
    -------
    Fqabc : np.ndarray of complex, shape (N,)
    """
    siA = sas_sinx_x(0.5*length_a*qa)
    siB = sas_sinx_x(0.5*length_b*qb)
    siC = sas_sinx_x(0.5*length_c*qc)
    return siA * siB * siC

def Iq(q, sld, sld_solvent, length_a, length_b, length_c):
    """
    Computes the scattering intensity I(q) of nanoprisms averaged over all orientations using the Fibonacci quadrature.
    The number of points on the sphere is set by npoints_fibonacci. Each point has an equal weight = 1/npoints_fibonacci.
    Parameters
    ----------
    q : float ou array
        Norm of the scattering vector
    sld, sld_solvent :
        Contrast of scattering length density
    Returns
    -------
    Iq : ndarray
        Scattering intensity averaged over all orientations
    """
    raise NotImplementedError()

Iq.vectorized = True

def Iqabc(qa, qb, qc, sld, sld_solvent, length_a, length_b, length_c):
    intensity = _fq(qa, qb, qc, length_a, length_b, length_c)
    volume = form_volume(length_a, length_b, length_c)
    #print("iqabc box", sld, sld_solvent, volume, length_a, length_b, length_c)
    s = 1e-2 * (sld - sld_solvent) * volume
    return (s * intensity)**2

tests = [
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_average": 10, "length": 200, "sld": 1., "sld_solvent": 0.},
     0.01, 5.62789],
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_average": 10, "length": 200, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [5.62789, 0.73696]],
]
