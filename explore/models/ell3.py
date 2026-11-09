"""
Pure python triaxial ellipsoid (2d only)
"""
import numpy as np
from numpy import inf

name = "ell3"
title = "pure python triaxial ellipsoid"
description = """
        Model for prisms of different cross-sections with orientation average using the Fibonacci quadrature"""
category = "shape:cylinder"
#             ["name", "units", default, [lower, upper], "type", "description"],
parameters = [["sld", "1e-6/Ang^2", 4, [-inf, inf], "sld",
               "Ellipsoid scattering length density"],
              ["sld_solvent", "1e-6/Ang^2", 1, [-inf, inf], "sld",
               "Solvent scattering length density"],
              ["radius_equat_minor", "Ang", 20, [0, inf], "volume",
               "Minor equatorial radius, Ra"],
              ["radius_equat_major", "Ang", 400, [0, inf], "volume",
               "Major equatorial radius, Rb"],
              ["radius_polar", "Ang", 10, [0, inf], "volume",
               "Polar radius, Rc"],
              ["theta", "degrees", 60, [-180, 180], "orientation",
               "polar axis to beam angle"],
              ["phi", "degrees", 60, [-360, 360], "orientation",
               "rotation about beam"],
              ["psi", "degrees", 60, [-360, 360], "orientation",
               "rotation about polar axis"],
             ]


### Functions for geometrical calculations: volume, surface of the cross section, average radius, edge length etc.
def form_volume(radius_equat_minor, radius_equat_major, radius_polar):
    """
    Computes the volume of a triaxial ellipsoid.

    Returns
    -------
    volume : float
    """
    return 4*np.pi/3*radius_equat_minor*radius_equat_major*radius_polar

def _fq(qa, qb, qc, radius_equat_minor, radius_equat_major, radius_polar):
    """
    Computes the form factor amplitude of a cylinder
    Parameters
    ----------
    qa, qb, qc : np.ndarray shape (N,)
        components of the scattering vector q
    Returns
    -------
    Fqabc : np.ndarray of real, shape (N,)
    """
    from sasmodels.special import sas_3j1x_x
    qr = np.sqrt(
        (radius_equat_minor*qa)**2
        + (radius_equat_major*qb)**2
        + (radius_polar*qc)**2
    )
    return sas_3j1x_x(qr)

def Iq(q, sld, sld_solvent, radius_equat_minor, radius_equat_major, radius_polar):
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

def Iqabc(qa, qb, qc, sld, sld_solvent, radius_equat_minor, radius_equat_major, radius_polar):
    intensity = _fq(qa, qb, qc, radius_equat_minor, radius_equat_major, radius_polar)
    volume = form_volume(radius_equat_minor, radius_equat_major, radius_polar)
    #print("iqabc box", sld, sld_solvent, volume, length_a, length_b, length_c)
    s = 1e-2 * (sld - sld_solvent) * volume
    return (s * intensity)**2

tests = [
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_average": 10, "length": 200, "sld": 1., "sld_solvent": 0.},
     0.01, 5.62789],
    [{"background": 0, "scale": 1, "n_sides": 4, "radius_average": 10, "length": 200, "sld": 1., "sld_solvent": 0.},
     [0.01, 0.1], [5.62789, 0.73696]],
]
