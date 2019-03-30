from __future__ import division, print_function
# Make sasmodels available on the path
import sys, os
BETA_DIR = os.path.dirname(os.path.realpath(__file__))
SASMODELS_DIR = os.path.dirname(os.path.dirname(BETA_DIR))
sys.path.insert(0, SASMODELS_DIR)

from collections import namedtuple

from matplotlib import pyplot as plt
import numpy as np
from numpy import pi, sin, cos, sqrt, fabs
from numpy.polynomial.legendre import leggauss
from scipy.special import j1 as J1
from numpy import inf
from scipy.special import gammaln  # type: ignore

Theory = namedtuple('Theory', 'Q F1 F2 P S I Seff Ibeta')
Theory.__new__.__defaults__ = (None,) * len(Theory._fields)

#Used to calculate F(q) for the cylinder, sphere, ellipsoid models
# RKH adding vesicle and hollow_cylinder to test sasview special cases of ER and VR
# There were issues here from python3 (i.e. qt5 version of sasview),fixed by Paul K (else do "source activate sasview")
def sas_sinx_x(x):
    with np.errstate(all='ignore'):
        retvalue = sin(x)/x
    retvalue[x == 0.] = 1.
    return retvalue

def sas_2J1x_x(x):
    with np.errstate(all='ignore'):
        retvalue = 2*J1(x)/x
    retvalue[x == 0] = 1.
    return retvalue

def sas_3j1x_x(x):
    """return 3*j1(x)/x"""
    retvalue = np.empty_like(x)
    with np.errstate(all='ignore'):
        # GSL bessel_j1 taylor expansion
        index = (x < 0.25)
        y = x[index]**2
        c1 = -1.0/10.0
        c2 = +1.0/280.0
        c3 = -1.0/15120.0
        c4 = +1.0/1330560.0
        c5 = -1.0/172972800.0
        retvalue[index] = 1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*c5))))
        index = ~index
        y = x[index]
        retvalue[index] = 3*(sin(y) - y*cos(y))/y**3
    retvalue[x == 0.] = 1.
    return retvalue

#Used to cross check my models with sasview models
def build_model(model_name, q, **pars):
    from sasmodels.core import load_model_info, build_model as build_sasmodel
    from sasmodels.data import empty_data1D
    from sasmodels.direct_model import DirectModel
    model_info = load_model_info(model_name)
    model = build_sasmodel(model_info, dtype='double!')
    data = empty_data1D(q)
    calculator = DirectModel(data, model,cutoff=0)
    calculator.pars = pars.copy()
    calculator.pars.setdefault('background', 0)
    return calculator

#gives the hardsphere structure factor that sasview uses
def _hardsphere_simple(q, radius_effective, volfraction):
    CUTOFFHS = 0.05
    if fabs(radius_effective) < 1.E-12:
        HARDSPH = 1.0
        return HARDSPH
    X = 1.0/(1.0 -volfraction)
    D = X*X
    A = (1.+2.*volfraction)*D
    A *= A
    X = fabs(q*radius_effective*2.0)
    if X < 5.E-06:
        HARDSPH = 1./A
        return HARDSPH
    X2 = X*X
    B = (1.0 +0.5*volfraction)*D
    B *= B
    B *= -6.*volfraction
    G = 0.5*volfraction*A
    if X < CUTOFFHS:
        FF = 8.0*A +6.0*B + 4.0*G + ( -0.8*A -B/1.5 -0.5*G +(A/35. +0.0125*B +0.02*G)*X2)*X2
        HARDSPH = 1./(1. + volfraction*FF )
        return HARDSPH
    X4 = X2*X2
    S, C = sin(X), cos(X)
    FF =  ((G*( (4.*X2 -24.)*X*S -(X4 -12.*X2 +24.)*C +24. )/X2 + B*(2.*X*S -(X2-2.)*C -2.) )/X + A*(S-X*C))/X
    HARDSPH = 1./(1. + 24.*volfraction*FF/X2)
    return HARDSPH

def hardsphere_simple(q, radius_effective, volfraction):
    SQ = [_hardsphere_simple(qk, radius_effective, volfraction) for qk in q]
    return np.array(SQ)

#Used in gaussian quadrature for polydispersity
#returns values and the probability of those values based on gaussian distribution
N_GAUSS = 35
NSIGMA_GAUSS = 3
def gaussian_distribution(center, sigma, lb, ub):
    #3 standard deviations covers approx. 99.7%
    if sigma != 0:
        nsigmas = NSIGMA_GAUSS
        x = np.linspace(center-sigma*nsigmas, center+sigma*nsigmas, num=N_GAUSS)
        x = x[(x >= lb) & (x <= ub)]
        px = np.exp((x-center)**2 / (-2.0 * sigma * sigma))
        return x, px
    else:
        return np.array([center]), np.array([1])

N_SCHULZ = 80
NSIGMA_SCHULZ = 8
def schulz_distribution(center, sigma, lb, ub):
    if sigma != 0:
        nsigmas = NSIGMA_SCHULZ
        x = np.linspace(center-sigma*nsigmas, center+sigma*nsigmas, num=N_SCHULZ)
        x = x[(x >= lb) & (x <= ub)]
        R = x/center
        z = (center/sigma)**2
        arg = z*np.log(z) + (z-1)*np.log(R) - R*z - np.log(center) - gammaln(z)
        px = np.exp(arg)
        return x, px
    else:
        return np.array([center]), np.array([1])

#returns the effective radius used in sasview
def ER_ellipsoid(radius_polar, radius_equatorial):
    ee = np.empty_like(radius_polar)
    if radius_polar > radius_equatorial:
        ee = (radius_polar**2 - radius_equatorial**2)/radius_polar**2
    elif radius_polar < radius_equatorial:
        ee = (radius_equatorial**2 - radius_polar**2) / radius_equatorial**2
    else:
        ee = 2*radius_polar
    if radius_polar * radius_equatorial != 0:
        bd = 1.0 - ee
        e1 = np.sqrt(ee)
        b1 = 1.0 + np.arcsin(e1) / (e1*np.sqrt(bd))
        bL = (1.0 + e1) / (1.0 - e1)
        b2 = 1.0 + bd / 2 / e1 * np.log(bL)
        delta = 0.75 * b1 * b2
    ddd = np.zeros_like(radius_polar)
    ddd = 2.0*(delta + 1.0)*radius_polar*radius_equatorial**2
    return 0.5*ddd**(1.0 / 3.0)

def ellipsoid_volume(radius_polar, radius_equatorial):
    volume = (4./3.)*pi*radius_polar*radius_equatorial**2
    return volume

# F1 is F(q)
# F2 is F(g)^2
#IQM is I(q) with monodispersity
#IQSM is I(q) with structure factor S(q) and monodispersity
#IQBM is I(q) with Beta Approximation and monodispersity
#SQ is monodisperse approach for structure factor
#SQ_EFF is the effective structure factor from beta approx
def ellipsoid_theta(q, radius_polar, radius_equatorial, sld, sld_solvent,
                    volfraction=0, radius_effective=None):
    #creates values z and corresponding probabilities w from legendre-gauss quadrature
    volume = ellipsoid_volume(radius_polar, radius_equatorial)
    z, w = leggauss(76)
    F1 = np.zeros_like(q)
    F2 = np.zeros_like(q)
    #use a u subsition(u=cos) and then u=(z+1)/2 to change integration from
    #0->2pi with respect to alpha to -1->1 with respect to z, allowing us to use
    #legendre-gauss quadrature
    for k, qk in enumerate(q):
        r = sqrt(radius_equatorial**2*(1-((z+1)/2)**2)+radius_polar**2*((z+1)/2)**2)
        form = (sld-sld_solvent)*volume*sas_3j1x_x(qk*r)
        F2[k] = np.sum(w*form**2)
        F1[k] = np.sum(w*form)
    #the 1/2 comes from the change of variables mentioned above
    F2 = F2/2.0
    F1 = F1/2.0
    if radius_effective is None:
        radius_effective = ER_ellipsoid(radius_polar,radius_equatorial)
    SQ = hardsphere_simple(q, radius_effective, volfraction)
    SQ_EFF = 1 + F1**2/F2*(SQ - 1)
    IQM = 1e-4*F2/volume
    IQSM = volfraction*IQM*SQ
    IQBM = volfraction*IQM*SQ_EFF
    return Theory(Q=q, F1=F1, F2=F2, P=IQM, S=SQ, I=IQSM, Seff=SQ_EFF, Ibeta=IQBM)

#IQD is I(q) polydispursed, IQSD is I(q)S(q) polydispursed, etc.
#IQBD HAS NOT BEEN CROSS CHECKED AT ALL
def ellipsoid_pe(q, radius_polar, radius_equatorial, sld, sld_solvent,
                 radius_polar_pd=0.1, radius_equatorial_pd=0.1,
                 radius_polar_pd_type='gaussian',
                 radius_equatorial_pd_type='gaussian',
                 volfraction=0, radius_effective=None,
                 background=0, scale=1,
                 norm='sasview'):
    if norm not in ['sasview', 'sasfit', 'yun']:
        raise TypeError("unknown norm "+norm)
    if radius_polar_pd_type == 'gaussian':
        Rp_val, Rp_prob = gaussian_distribution(radius_polar, radius_polar_pd*radius_polar, 0, inf)
    elif radius_polar_pd_type == 'schulz':
        Rp_val, Rp_prob = schulz_distribution(radius_polar, radius_polar_pd*radius_polar, 0, inf)
    if radius_equatorial_pd_type == 'gaussian':
        Re_val, Re_prob = gaussian_distribution(radius_equatorial, radius_equatorial_pd*radius_equatorial, 0, inf)
    elif radius_equatorial_pd_type == 'schulz':
        Re_val, Re_prob = schulz_distribution(radius_equatorial, radius_equatorial_pd*radius_equatorial, 0, inf)
    total_weight = total_volume = 0
    radius_eff = 0
    F1, F2 = np.zeros_like(q), np.zeros_like(q)
    for k, Rpk in enumerate(Rp_val):
        print("ellipsoid cycle", k, "of", len(Rp_val))
        for i, Rei in enumerate(Re_val):
            theory = ellipsoid_theta(q, Rpk, Rei, sld, sld_solvent)
            volume = ellipsoid_volume(Rpk, Rei)
            weight = Rp_prob[k]*Re_prob[i]
            total_weight += weight
            total_volume += weight*volume
            F1 += theory.F1*weight
            F2 += theory.F2*weight
            radius_eff += weight*ER_ellipsoid(Rpk, Rei)
    F1 /= total_weight
    F2 /= total_weight
    average_volume = total_volume/total_weight
    if radius_effective is None:
        radius_effective = radius_eff/total_weight
    if norm == 'sasfit':
        IQD = F2
    elif norm == 'sasview':
        # Note: internally, sasview uses F2/total_volume because:
        #   average_volume = total_volume/total_weight
        #   F2/total_weight / average_volume
        #     = F2/total_weight / total_volume/total_weight
        #     = F2/total_volume
        IQD = F2/average_volume*1e-4*volfraction
        F1 *= 1e-2  # Yun is using sld in 1/A^2, not 1e-6/A^2
        F2 *= 1e-4
    elif norm == 'yun':
        F1 *= 1e-6  # Yun is using sld in 1/A^2, not 1e-6/A^2
        F2 *= 1e-12
        IQD = F2/average_volume*1e8*volfraction
    SQ = hardsphere_simple(q, radius_effective, volfraction)
    beta = F1**2/F2
    SQ_EFF = 1 + beta*(SQ - 1)
    IQSD = IQD*SQ
    IQBD = IQD*SQ_EFF
    return Theory(Q=q, F1=F1, F2=F2, P=IQD, S=SQ, I=IQSD, Seff=SQ_EFF, Ibeta=IQBD)

#polydispersity for sphere
def sphere_r(q,radius,sld,sld_solvent,
             radius_pd=0.1, radius_pd_type='gaussian',
             volfraction=0, radius_effective=None,
             background=0, scale=1,
             norm='sasview'):
    if norm not in ['sasview', 'sasfit', 'yun']:
        raise TypeError("unknown norm "+norm)
    if radius_pd_type == 'gaussian':
        radius_val, radius_prob = gaussian_distribution(radius, radius_pd*radius, 0, inf)
    elif radius_pd_type == 'schulz':
        radius_val, radius_prob = schulz_distribution(radius, radius_pd*radius, 0, inf)
    total_weight = total_volume = 0
    F1 = np.zeros_like(q)
    F2 = np.zeros_like(q)
    for k, rk in enumerate(radius_val):
        volume = 4./3.*pi*rk**3
        total_weight += radius_prob[k]
        total_volume += radius_prob[k]*volume
        form = (sld-sld_solvent)*volume*sas_3j1x_x(q*rk)
        F2 += radius_prob[k]*form**2
        F1 += radius_prob[k]*form
    F1 /= total_weight
    F2 /= total_weight
    average_volume = total_volume/total_weight

    if radius_effective is None:
        radius_effective = radius
        average_volume = total_volume/total_weight
    if norm == 'sasfit':
        IQD = F2
    elif norm == 'sasview':
        IQD = F2/average_volume*1e-4*volfraction
    elif norm == 'yun':
        F1 *= 1e-6  # Yun is using sld in 1/A^2, not 1e-6/A^2
        F2 *= 1e-12
        IQD = F2/average_volume*1e8*volfraction
    print("in sphere_r : radius_effective  ",radius_effective,"  volfraction",volfraction)
    SQ = hardsphere_simple(q, radius_effective, volfraction)
    beta = F1**2/F2
    SQ_EFF = 1 + beta*(SQ - 1)
    IQSD = IQD*SQ
    IQBD = IQD*SQ_EFF
    return Theory(Q=q, F1=F1, F2=F2, P=IQD, S=SQ, I=IQSD, Seff=SQ_EFF, Ibeta=IQBD)

#polydispersity for vesicle
def vesicle_pe(q,radius=20, thickness=10, sld=4, sld_solvent=1, volfraction=0.03,
        radius_pd=0.1, thickness_pd=0.2, radius_pd_type="gaussian", thickness_pd_type="gaussian",
        radius_effective=None, background=0, scale=1,
                 norm='sasview'):
    if norm not in ['sasview', 'sasfit', 'yun']:
        raise TypeError("unknown norm "+norm)
    if radius_pd_type == 'gaussian':
        R_val, R_prob = gaussian_distribution(radius, radius_pd*radius, 0, inf)
    elif radius_pd_type == 'schulz':
        R_val, R_prob = schulz_distribution(radius, radius_pd*radius, 0, inf)
    if thickness_pd_type == 'gaussian':
        T_val, T_prob = gaussian_distribution(thickness, thickness_pd*thickness, 0, inf)
    elif thickness_pd_type == 'schulz':
        T_val, T_prob = schulz_distribution(thickness, thickness_pd*thickness, 0, inf)
    total_weight = total_volume = total_shell = 0
    radius_eff = 0
    F1, F2 = np.zeros_like(q), np.zeros_like(q)
    for k, Rk in enumerate(R_val):
        print("vesicle cycle", k, "of", len(R_val)," Rk ",Rk)
        volume_in = 4./3.*pi*Rk**3
        form_in = (sld-sld_solvent)*volume_in*sas_3j1x_x(q*Rk)
        for i, Ti in enumerate(T_val):
            Rout = Rk + Ti
            volume_out = 4./3.*pi*Rout**3
            form_out = (sld-sld_solvent)*volume_out*sas_3j1x_x(q*Rout)
            form = form_out - form_in
            volume = volume_out - volume_in
            weight = R_prob[k]*T_prob[i]
            total_weight += weight
            total_shell += weight*volume
            total_volume += weight*volume_out
            F1 += form*weight
            F2 += form**2*weight
            radius_eff += weight*Rout
    F1 /= total_weight
    F2 /= total_weight
    average_volume = total_volume/total_weight
    if radius_effective is None:
        radius_effective = radius_eff/total_weight
        print("in vesicle_pe : radius_effective  ",radius_effective)
    if norm == 'sasfit':
        IQD = F2
    elif norm == 'sasview':
        # Note: internally, sasview uses F2/total_volume because:
        #   average_volume = total_volume/total_weight
        #   F2/total_weight / average_volume
        #     = F2/total_weight / total_volume/total_weight
        #     = F2/total_volume
        IQD = F2/average_volume*1e-4*volfraction
        F1 *= 1e-2  # Yun is using sld in 1/A^2, not 1e-6/A^2
        F2 *= 1e-4
    elif norm == 'yun':
        F1 *= 1e-6  # Yun is using sld in 1/A^2, not 1e-6/A^2
        F2 *= 1e-12
        IQD = F2/average_volume*1e8*volfraction
    # RKH SORT THIS OUT WHEN HAVE NEW MODEL  Vtot/Vshell = (R+T)^3/ ( (R+T)^3 -R^3 )
    vfff=volfraction*( 1.0/(1.0 - ( (radius/(radius+thickness))**3 )))
    print("in vesicle_pe : radius_effective, fudged volfraction for S(Q) ",radius_effective,vfff)
    SQ = hardsphere_simple(q, radius_effective, vfff)
    beta = F1**2/F2
    SQ_EFF = 1 + beta*(SQ - 1)
    IQSD = IQD*SQ
    IQBD = IQD*SQ_EFF
    return Theory(Q=q, F1=F1, F2=F2, P=IQD, S=SQ, I=IQSD, Seff=SQ_EFF, Ibeta=IQBD)
#
#polydispersity for hollow_cylinder
def hollow_cylinder_pe(q,radius=20, thickness=10, length=80, sld=4, sld_solvent=1, volfraction=0.15,
        radius_pd=0.1, thickness_pd=0.2, length_pd=0.05, radius_pd_type="gaussian", length_pd_type="gaussian",
        thickness_pd_type="gaussian", radius_effective=None, background=0, scale=1,
                 norm='sasview'):
    if norm not in ['sasview', 'sasfit', 'yun']:
        raise TypeError("unknown norm "+norm)
    if radius_pd_type == 'gaussian':
        R_val, R_prob = gaussian_distribution(radius, radius_pd*radius, 0, inf)
    elif radius_pd_type == 'schulz':
        R_val, R_prob = schulz_distribution(radius, radius_pd*radius, 0, inf)
    if thickness_pd_type == 'gaussian':
        T_val, T_prob = gaussian_distribution(thickness, thickness_pd*thickness, 0, inf)
    elif thickness_pd_type == 'schulz':
        T_val, T_prob = schulz_distribution(thickness, thickness_pd*thickness, 0, inf)
    if length_pd_type == 'gaussian':
        L_val, L_prob = gaussian_distribution(length, length_pd*length, 0, inf)
    elif length_pd_type == 'schulz':
        L_val, L_prob = schulz_distribution(length, length_pd*length, 0, inf)
    total_weight = total_volume = total_shell = 0
    radius_eff = 0
    F1, F2 = np.zeros_like(q), np.zeros_like(q)
    for k, Rk in enumerate(R_val):
        print("hollow_cylinder cycle", k, "of", len(R_val)," Rk ",Rk)
        for i, Ti in enumerate(T_val):
            for j, Lj in enumerate(L_val):
                Rout = Rk + Ti
                volume_out = pi*Rout**2*Lj
                volume_in = pi*Rk**2*Lj
                volume = volume_out - volume_in
                theory = hollow_cylinder_theta(q, Rk, Ti, Lj, sld, sld_solvent,
                    volfraction=0, radius_effective=None)
                weight = R_prob[k]*T_prob[i]*L_prob[j]
                total_weight += weight
                total_shell += weight*volume
                total_volume += weight*volume_out
                F1 += theory.F1*weight
                F2 += theory.F2*weight
                radius_eff += weight*ER_hollow_cylinder(radius, thickness, length)
    F1 /= total_weight
    F2 /= total_weight
    average_volume = total_volume/total_weight
    if radius_effective is None:
        radius_effective = radius_eff/total_weight
        print("in hollow_cylinder : radius_effective  ",radius_effective)
    if norm == 'sasfit':
        IQD = F2
    elif norm == 'sasview':
        # Note: internally, sasview uses F2/total_volume because:
        #   average_volume = total_volume/total_weight
        #   F2/total_weight / average_volume
        #     = F2/total_weight / total_volume/total_weight
        #     = F2/total_volume
        IQD = F2/average_volume*1e-4*volfraction
        F1 *= 1e-2  # Yun is using sld in 1/A^2, not 1e-6/A^2
        F2 *= 1e-4
    elif norm == 'yun':
        F1 *= 1e-6  # Yun is using sld in 1/A^2, not 1e-6/A^2
        F2 *= 1e-12
        IQD = F2/average_volume*1e8*volfraction
# RKH SORT THIS OUT WHEN HAVE NEW MODEL
    vfff=volfraction
    print("in hollow_cylinder : radius_effective, volfaction for S(Q) ",radius_effective,vfff)
    SQ = hardsphere_simple(q, radius_effective, vfff)
    beta = F1**2/F2
    SQ_EFF = 1 + beta*(SQ - 1)
    IQSD = IQD*SQ
    IQBD = IQD*SQ_EFF
    return Theory(Q=q, F1=F1, F2=F2, P=IQD, S=SQ, I=IQSD, Seff=SQ_EFF, Ibeta=IQBD)

# RKH copying Paul K's ellipsoid here, this computes everything for monodisperse case
# and supplies partial results for polydisperse case.
def hollow_cylinder_theta(q, radius, thickness, length, sld, sld_solvent,
                    volfraction=0, radius_effective=None):
    #creates values z and corresponding probabilities w from legendre-gauss quadrature
    z, w = leggauss(76)
    F1 = np.zeros_like(q)
    F2 = np.zeros_like(q)
    #use a u subsition(u=cos) and then u=(z+1)/2 to change integration from
    #0->2pi with respect to alpha to -1->1 with respect to z, allowing us to use
    #legendre-gauss quadrature
    lower = 0.0
    upper = 1.0
    gamma_sq = (radius/(radius+thickness))**2
    for k, qk in enumerate(q):
        for i, zt in enumerate(z):
            # quadrature loop
            cos_theta = 0.5*( zt* (upper-lower) + lower + upper  )
            sin_theta = sqrt(1.0 - cos_theta*cos_theta)
            aaa=(radius+thickness)*qk*sin_theta
            ## sas_J1 expects an array, so try direct use of J1
            lam1 = J1(aaa)/aaa
            aaa=radius*qk*sin_theta
            lam2 = J1(aaa)/aaa
            psi = (lam1 - gamma_sq*lam2)/(1.0 - gamma_sq)
            #Note: lim_{thickness -> 0} psi = sas_J0(radius*qab)
            #Note: lim_{radius -> 0} psi = sas_2J1x_x(thickness*qab)
            aaa=0.5*length*qk*cos_theta
            t2 = sin(aaa)/aaa
#            t2 = sas_sinx_x(0.5*length*qk*cos_theta)
            form = psi*t2
            F1[i] += w[i] * form
            F2[i] += w[i] * form * form
    volume = pi*((radius + thickness)**2 - radius**2)*length
    s = (sld - sld_solvent) * volume
    F2 = F2*s*(upper - lower)/2.0
    F1 = F1*s*s*(upper - lower)/2.0
    if radius_effective is None:
        radius_effective = ER_hollow_cylinder(radius, thickness, length)
    SQ = hardsphere_simple(q, radius_effective, volfraction)
    SQ_EFF = 1 + F1**2/F2*(SQ - 1)
    IQM = 1e-4*F2/volume
    IQSM = volfraction*IQM*SQ
    IQBM = volfraction*IQM*SQ_EFF
    return Theory(Q=q, F1=F1, F2=F2, P=IQM, S=SQ, I=IQSM, Seff=SQ_EFF, Ibeta=IQBM)
#
def ER_hollow_cylinder(radius, thickness, length):
    """
    :param radius:      Cylinder core radius
    :param thickness:   Cylinder wall thickness
    :param length:      Cylinder length
    :return:            Effective radius
    """
    router = radius + thickness
    if router == 0 or length == 0:
        return 0.0
    len1 = router
    len2 = length/2.0
    term1 = len1*len1*2.0*len2/2.0
    term2 = 1.0 + (len2/len1)*(1.0 + 1/len2/2.0)*(1.0 + pi*len1/len2/2.0)
    ddd = 3.0*term1*term2
    diam = pow(ddd, (1.0/3.0))
    return diam

###############################################################################
###############################################################################
###############################################################################
##################                                           ##################
##################                   TESTS                   ##################
##################                                           ##################
###############################################################################
###############################################################################
###############################################################################

def popn(d, keys):
    """
    Splits a dict into two, with any key of *d* which is in *keys* removed
    from *d* and added to *b*. Returns *b*.
    """
    b = {}
    for k in keys:
        try:
            b[k] = d.pop(k)
        except KeyError:
            pass
    return b

def sasmodels_theory(q, Pname, **pars):
    """
    Call sasmodels to compute the model with and without a hard sphere
    structure factor.
    """
    #mono_pars = {k: (0 if k.endswith('_pd') else v) for k, v in pars.items()}
    Ppars = pars.copy()
    Spars = popn(Ppars, ['radius_effective', 'volfraction'])
    Ipars = pars.copy()

    # Autofill npts and nsigmas for the given pd type
    for k, v in pars.items():
        if k.endswith("_pd_type"):
            if v == "gaussian":
                n, nsigmas = N_GAUSS, NSIGMA_GAUSS
            elif v == "schulz":
                n, nsigmas = N_SCHULZ, NSIGMA_SCHULZ
            Ppars.setdefault(k.replace("_pd_type", "_pd_n"), n)
            Ppars.setdefault(k.replace("_pd_type", "_pd_nsigma"), nsigmas)
            Ipars.setdefault(k.replace("_pd_type", "_pd_n"), n)
            Ipars.setdefault(k.replace("_pd_type", "_pd_nsigma"), nsigmas)

    #Ppars['scale'] = Spars.get('volfraction', 1)
    P = build_model(Pname, q)
    S = build_model("hardsphere", q)
    I = build_model(Pname+"@hardsphere", q)
    Pq = P(**Ppars)*pars.get('volfraction', 1)
    Sq = S(**Spars)
    Iq = I(**Ipars)
    #Iq = Pq*Sq*pars.get('volfraction', 1)
    #Sq = Iq/Pq
    #Iq = None#= Sq = None
    r = dict(I._kernel.results())
    print(r)
    return Theory(Q=q, F1=None, F2=None, P=Pq, S=None, I=None, Seff=r["S_eff(Q)"], Ibeta=Iq)

def compare(title, target, actual, fields='F1 F2 P S I Seff Ibeta'):
    """
    Plot fields in common between target and actual, along with relative error.
    """
    available = [s for s in fields.split()
                 if getattr(target, s) is not None and getattr(actual, s) is not None]
    rows = len(available)
    for row, field in enumerate(available):
        Q = target.Q
        I1, I2 = getattr(target, field), getattr(actual, field)
        plt.subplot(rows, 2, 2*row+1)
        plt.loglog(Q, abs(I1), label="target "+field)
        plt.loglog(Q, abs(I2), label="value "+field)
        #plt.legend(loc="upper left", bbox_to_anchor=(1,1))
        plt.legend(loc='lower left')
        plt.subplot(rows, 2, 2*row+2)
        plt.semilogx(Q, I2/I1 - 1, label="relative error")
        #plt.semilogx(Q, I1/I2 - 1, label="relative error")
    plt.tight_layout()
    plt.suptitle(title)
    plt.show()

def data_file(name):
    return os.path.join(BETA_DIR, 'data_files', name)

def load_sasfit(path):
    data = np.loadtxt(path, dtype=str, delimiter=';').T
    data = np.vstack((list(map(float, v)) for v in data[0:2]))
    return data

COMPARISON = {}  # Type: Dict[(str,str,str)] -> Callable[(), None]

def compare_sasview_sphere(pd_type='schulz'):
    q = np.logspace(-5, 0, 250)
    model = 'sphere'
    pars = dict(
        radius=20, sld=4, sld_solvent=1,
        background=0,
        radius_pd=.1, radius_pd_type=pd_type,
        volfraction=0.15,
        radius_effective=20  # equivalent average sphere radius
#       NOTE sasview computes its own radius_effective in "target" (the print(r) at end of sasmodels_theory will reveal its value),
#       this one is only used locally for "actual"
        )
    target = sasmodels_theory(q, model, radius_effective_mode=0, structure_factor_mode=1, **pars)
    actual = sphere_r(q, norm='sasview', **pars)
    title = " ".join(("sasmodels", model, pd_type))
    compare(title, target, actual)
COMPARISON[('sasview', 'sphere', 'gaussian')] = lambda: compare_sasview_sphere(pd_type='gaussian')
COMPARISON[('sasview', 'sphere', 'schulz')] = lambda: compare_sasview_sphere(pd_type='schulz')


def compare_sasview_vesicle(pd_type='gaussian'):
    q = np.logspace(-5, 0, 250)
    model = 'vesicle'
    print("F2 seems OK, but big issues with S(Q), so work in progress")
# NOTE parameters for vesicle model are soon to be changed to remove "volfraction" (though still there in 5.0)
    pars = dict(
        radius=20, thickness=10, sld=4, sld_solvent=1, volfraction=0.03,
        background=0,
        radius_pd=0.01, thickness_pd=0.01, radius_pd_type=pd_type, thickness_pd_type=pd_type,
        radius_effective=30. )
        # equivalent average sphere radius for local "actual" to match what sasview uses, use None to compute average outer radius here,

    target = sasmodels_theory(q, model, radius_effective_mode=0, structure_factor_mode=1, **pars)
    actual = vesicle_pe(q, norm='sasview', **pars)
    title = " ".join(("sasmodels", model, pd_type))
    compare(title, target, actual)
COMPARISON[('sasview', 'vesicle', 'gaussian')] = lambda: compare_sasview_vesicle(pd_type='gaussian')
COMPARISON[('sasview', 'vesicle', 'schulz')] = lambda: compare_sasview_vesicle(pd_type='schulz')

def compare_sasview_hollow_cylinder(pd_type='gaussian'):
    q = np.logspace(-5, 0, 250)
    model = 'hollow_cylinder'
    print(model)
    print("just about works for monodisperse, polydisperse needs work")
# NOTE parameters for hollow_cylinder model are soon to be changed to remove "volfraction"
# setting all three pd to zero is OK
    pars = dict(
        radius=20, thickness=10, length=80, sld=4, sld_solvent=1,
        background=0,
        radius_pd=0.1, thickness_pd=0.0, length_pd=0.0, radius_pd_type=pd_type, thickness_pd_type=pd_type, length_pd_type=pd_type,
        radius_effective=40.687)
        # equivalent average sphere radius for local "actual" to match what sasview uses
    target = sasmodels_theory(q, model, radius_effective_mode=0, structure_factor_mode=1, **pars)
    actual = hollow_cylinder_pe(q, norm='sasview', **pars)
# RKH monodisp was OK,    actual = hollow_cylinder_theta(q,radius=20, thickness=10, length=80, sld=4, sld_solvent=1 )
    title = " ".join(("sasmodels", model, pd_type))
    compare(title, target, actual)
COMPARISON[('sasview', 'hollow_cylinder', 'gaussian')] = lambda: compare_sasview_hollow_cylinder(pd_type='gaussian')
COMPARISON[('sasview', 'hollow_cylinder', 'schulz')] = lambda: compare_sasview_hollow_cylinder(pd_type='schulz')


def compare_sasview_ellipsoid(pd_type='gaussian'):
    q = np.logspace(-5, 0, 50)
    model = 'ellipsoid'
    pars = dict(
        radius_polar=20, radius_equatorial=400, sld=4, sld_solvent=1,
        background=0,
        radius_polar_pd=0.1, radius_polar_pd_type=pd_type,
        radius_equatorial_pd=0.1, radius_equatorial_pd_type=pd_type,
        volfraction=0.15,
        radius_effective=270.7543927018,
# if change radius_effective to some other value, the S(Q) from sasview does not agree
        )
    target = sasmodels_theory(q, model, radius_effective_mode=0, structure_factor_mode=1, **pars)
    actual = ellipsoid_pe(q, norm='sasview', **pars)
# RKH test       actual = ellipsoid_theta(q, radius_polar=20, radius_equatorial=400, sld=4, sld_solvent=1, volfraction=0.15, radius_effective=270.)
    title = " ".join(("sasmodels", model, pd_type))
    compare(title, target, actual)
COMPARISON[('sasview', 'ellipsoid', 'gaussian')] = lambda: compare_sasview_ellipsoid(pd_type='gaussian')
COMPARISON[('sasview', 'ellipsoid', 'schulz')] = lambda: compare_sasview_ellipsoid(pd_type='schulz')

def compare_yun_ellipsoid_mono():
    pars = {
        'radius_polar': 20, 'radius_polar_pd': 0, 'radius_polar_pd_type': 'gaussian',
        'radius_equatorial': 10, 'radius_equatorial_pd': 0, 'radius_equatorial_pd_type': 'gaussian',
        'sld': 2, 'sld_solvent': 1,
        'volfraction': 0.15,
        # Yun uses radius for same volume sphere for effective radius
        # whereas sasview uses the average curvature.
        'radius_effective': 12.59921049894873,
    }

    data = np.loadtxt(data_file('yun_ellipsoid.dat'),skiprows=2).T
    Q = data[0]
    F1 = data[1]
    P = data[3]*pars['volfraction']
    S = data[5]
    Seff = data[6]
    target = Theory(Q=Q, F1=F1, P=P, S=S, Seff=Seff)
    actual = ellipsoid_pe(Q, norm='yun', **pars)
    title = " ".join(("yun", "ellipsoid", "no pd"))
    #compare(title, target, actual, fields="P S I Seff Ibeta")
    compare(title, target, actual)
COMPARISON[('yun', 'ellipsoid', 'gaussian')] = compare_yun_ellipsoid_mono
COMPARISON[('yun', 'ellipsoid', 'schulz')] = compare_yun_ellipsoid_mono

def compare_yun_sphere_gauss():
    # Note: yun uses gauss limits from R0/10 to R0 + 5 sigma steps sigma/100
    # With pd = 0.1, that's 14 sigma and 1400 points.
    pars = {
        'radius': 20, 'radius_pd': 0.1, 'radius_pd_type': 'gaussian',
        'sld': 6, 'sld_solvent': 0,
        'volfraction': 0.1,
    }

    data = np.loadtxt(data_file('testPolydisperseGaussianSphere.dat'), skiprows=2).T
    Q = data[0]
    F1 = data[1]
    F2 = data[2]
    P = data[3]
    S = data[5]
    Seff = data[6]
    target = Theory(Q=Q, F1=F1, P=P, S=S, Seff=Seff)
    actual = sphere_r(Q, norm='yun', **pars)
    title = " ".join(("yun", "sphere", "10% dispersion 10% Vf"))
    compare(title, target, actual)
    data = np.loadtxt(data_file('testPolydisperseGaussianSphere2.dat'), skiprows=2).T
    pars.update(radius_pd=0.15)
    Q = data[0]
    F1 = data[1]
    F2 = data[2]
    P = data[3]
    S = data[5]
    Seff = data[6]
    target = Theory(Q=Q, F1=F1, P=P, S=S, Seff=Seff)
    actual = sphere_r(Q, norm='yun', **pars)
    title = " ".join(("yun", "sphere", "15% dispersion 10% Vf"))
    compare(title, target, actual)
COMPARISON[('yun', 'sphere', 'gaussian')] = compare_yun_sphere_gauss


def compare_sasfit_sphere_gauss():
    #N=1,s=2,X0=20,distr radius R=20,eta_core=4,eta_solv=1,.3
    pars = {
        'radius': 20, 'radius_pd': 0.1, 'radius_pd_type': 'gaussian',
        'sld': 4, 'sld_solvent': 1,
        'volfraction': 0.3,
    }

    Q, IQ = load_sasfit(data_file('sasfit_sphere_IQD.txt'))
    Q, IQSD = load_sasfit(data_file('sasfit_sphere_IQSD.txt'))
    Q, IQBD = load_sasfit(data_file('sasfit_sphere_IQBD.txt'))
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_sphere_sq.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_sphere_sqeff.txt'))
    target = Theory(Q=Q, F1=None, F2=None, P=IQ, S=SQ, I=IQSD, Seff=SQ_EFF, Ibeta=IQBD)
    actual = sphere_r(Q, norm="sasfit", **pars)
    title = " ".join(("sasfit", "sphere", "pd=10% gaussian"))
    compare(title, target, actual)
    #compare(title, target, actual, fields="P")
COMPARISON[('sasfit', 'sphere', 'gaussian')] = compare_sasfit_sphere_gauss

def compare_sasfit_sphere_schulz():
    #radius=20,sld=4,sld_solvent=1,volfraction=0.3,radius_pd=0.1
    #We have scaled the output from sasfit by 1e-4*volume*volfraction
    #0.10050378152592121
    pars = {
        'radius': 20, 'radius_pd': 0.1, 'radius_pd_type': 'schulz',
        'sld': 4, 'sld_solvent': 1,
        'volfraction': 0.3,
    }

    Q, IQ = load_sasfit(data_file('sasfit_sphere_schulz_IQD.txt'))
    Q, IQSD = load_sasfit(data_file('sasfit_sphere_schulz_IQSD.txt'))
    Q, IQBD = load_sasfit(data_file('sasfit_sphere_schulz_IQBD.txt'))
    target = Theory(Q=Q, F1=None, F2=None, P=IQ, S=None, I=IQSD, Seff=None, Ibeta=IQBD)
    actual = sphere_r(Q, norm="sasfit", **pars)
    title = " ".join(("sasfit", "sphere", "pd=10% schulz"))
    compare(title, target, actual)
COMPARISON[('sasfit', 'sphere', 'schulz')] = compare_sasfit_sphere_schulz

def compare_sasfit_ellipsoid_schulz():
    #polarradius=20, equatorialradius=10, sld=4,sld_solvent=1,volfraction=0.3,radius_polar_pd=0.1
    #Effective radius =13.1353356684
    #We have scaled the output from sasfit by 1e-4*volume*volfraction
    #0.10050378152592121
    pars = {
        'radius_polar': 20, 'radius_polar_pd': 0.1, 'radius_polar_pd_type': 'schulz',
        'radius_equatorial': 10, 'radius_equatorial_pd': 0., 'radius_equatorial_pd_type': 'schulz',
        'sld': 4, 'sld_solvent': 1,
        'volfraction': 0.3, 'radius_effective': 13.1353356684,
    }

    Q, IQ = load_sasfit(data_file('sasfit_ellipsoid_schulz_IQD.txt'))
    Q, IQSD = load_sasfit(data_file('sasfit_ellipsoid_schulz_IQSD.txt'))
    Q, IQBD = load_sasfit(data_file('sasfit_ellipsoid_schulz_IQBD.txt'))
    target = Theory(Q=Q, F1=None, F2=None, P=IQ, S=None, I=IQSD, Seff=None, Ibeta=IQBD)
    actual = ellipsoid_pe(Q, norm="sasfit", **pars)
    title = " ".join(("sasfit", "ellipsoid", "pd=10% schulz"))
    compare(title, target, actual)
COMPARISON[('sasfit', 'ellipsoid', 'schulz')] = compare_sasfit_ellipsoid_schulz


def compare_sasfit_ellipsoid_gaussian():
    pars = {
        'radius_polar': 20, 'radius_polar_pd': 0, 'radius_polar_pd_type': 'gaussian',
        'radius_equatorial': 10, 'radius_equatorial_pd': 0, 'radius_equatorial_pd_type': 'gaussian',
        'sld': 4, 'sld_solvent': 1,
        'volfraction': 0, 'radius_effective': None,
    }

    #Rp=20,Re=10,eta_core=4,eta_solv=1
    Q, PQ0 = load_sasfit(data_file('sasfit_ellipsoid_IQM.txt'))
    pars.update(volfraction=0, radius_polar_pd=0.0, radius_equatorial_pd=0, radius_effective=None)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, P=PQ0)
    compare("sasfit ellipsoid no poly", target, actual); plt.show()

    #N=1,s=2,X0=20,distr 10% polar Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    Q, PQ_Rp10 = load_sasfit(data_file('sasfit_ellipsoid_IQD.txt'))
    pars.update(volfraction=0, radius_polar_pd=0.1, radius_equatorial_pd=0.0, radius_effective=None)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, P=PQ_Rp10)
    compare("sasfit ellipsoid P(Q) 10% Rp", target, actual); plt.show()
    #N=1,s=1,X0=10,distr 10% equatorial Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    Q, PQ_Re10 = load_sasfit(data_file('sasfit_ellipsoid_IQD2.txt'))
    pars.update(volfraction=0, radius_polar_pd=0.0, radius_equatorial_pd=0.1, radius_effective=None)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, P=PQ_Re10)
    compare("sasfit ellipsoid P(Q) 10% Re", target, actual); plt.show()
    #N=1,s=6,X0=20,distr 30% polar Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    Q, PQ_Rp30 = load_sasfit(data_file('sasfit_ellipsoid_IQD3.txt'))
    pars.update(volfraction=0, radius_polar_pd=0.3, radius_equatorial_pd=0.0, radius_effective=None)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, P=PQ_Rp30)
    compare("sasfit ellipsoid P(Q) 30% Rp", target, actual); plt.show()
    #N=1,s=3,X0=10,distr 30% equatorial Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    Q, PQ_Re30 = load_sasfit(data_file('sasfit_ellipsoid_IQD4.txt'))
    pars.update(volfraction=0, radius_polar_pd=0.0, radius_equatorial_pd=0.3, radius_effective=None)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, P=PQ_Re30)
    compare("sasfit ellipsoid P(Q) 30% Re", target, actual); plt.show()
    #N=1,s=12,X0=20,distr 60% polar Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    Q, PQ_Rp60 = load_sasfit(data_file('sasfit_ellipsoid_IQD5.txt'))
    pars.update(volfraction=0, radius_polar_pd=0.6, radius_equatorial_pd=0.0, radius_effective=None)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, P=PQ_Rp60)
    compare("sasfit ellipsoid P(Q) 60% Rp", target, actual); plt.show()
    #N=1,s=6,X0=10,distr 60% equatorial Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    Q, PQ_Re60 = load_sasfit(data_file('sasfit_ellipsoid_IQD6.txt'))
    pars.update(volfraction=0, radius_polar_pd=0.0, radius_equatorial_pd=0.6, radius_effective=None)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, P=PQ_Re60)
    compare("sasfit ellipsoid P(Q) 60% Re", target, actual); plt.show()

    #N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.15
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff.txt'))
    pars.update(volfraction=0.15, radius_polar_pd=0.1, radius_equatorial_pd=0, radius_effective=13.1354236254)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 10% Rp 15% Vf", target, actual); plt.show()
    #N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.15
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq2.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff2.txt'))
    pars.update(volfraction=0.15, radius_polar_pd=0.3, radius_equatorial_pd=0, radius_effective=13.0901197149)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 30% Rp 15% Vf", target, actual); plt.show()
    #N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.15
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq3.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff3.txt'))
    pars.update(volfraction=0.15, radius_polar_pd=0.6, radius_equatorial_pd=0, radius_effective=13.336060917)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 60% Rp 15% Vf", target, actual); plt.show()

    #N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.3
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq4.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff4.txt'))
    pars.update(volfraction=0.3, radius_polar_pd=0.1, radius_equatorial_pd=0, radius_effective=13.1354236254)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 10% Rp 30% Vf", target, actual); plt.show()
    #N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.3
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq5.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff5.txt'))
    pars.update(volfraction=0.3, radius_polar_pd=0.3, radius_equatorial_pd=0, radius_effective=13.0901197149)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 30% Rp 30% Vf", target, actual); plt.show()
    #N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.3
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq6.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff6.txt'))
    pars.update(volfraction=0.3, radius_polar_pd=0.6, radius_equatorial_pd=0, radius_effective=13.336060917)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 60% Rp 30% Vf", target, actual); plt.show()

    #N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.6
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq7.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff7.txt'))
    pars.update(volfraction=0.6, radius_polar_pd=0.1, radius_equatorial_pd=0, radius_effective=13.1354236254)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 10% Rp 60% Vf", target, actual); plt.show()
    #N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.6
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq8.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff8.txt'))
    pars.update(volfraction=0.6, radius_polar_pd=0.3, radius_equatorial_pd=0, radius_effective=13.0901197149)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 30% Rp 60% Vf", target, actual); plt.show()
    #N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.6
    Q, SQ = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sq9.txt'))
    Q, SQ_EFF = load_sasfit(data_file('sasfit_polydisperse_ellipsoid_sqeff9.txt'))
    pars.update(volfraction=0.6, radius_polar_pd=0.6, radius_equatorial_pd=0, radius_effective=13.336060917)
    actual = ellipsoid_pe(Q, norm='sasfit', **pars)
    target = Theory(Q=Q, S=SQ, Seff=SQ_EFF)
    compare("sasfit ellipsoid P(Q) 60% Rp 60% Vf", target, actual); plt.show()
COMPARISON[('sasfit', 'ellipsoid', 'gaussian')] = compare_sasfit_ellipsoid_gaussian

def main():
    key = tuple(sys.argv[1:])
    if key not in COMPARISON:
        print("Usage: python sasfit_compare.py [sasview|sasfit|yun] [sphere|ellipsoid|vesicle|hollow_cylinder] [gaussian|schulz]")
        print("But not for [yun sphere schulz], and vesicle & hollow_cylinder only with sasview.")
        print("Note this compares chosen source against internal python code here, with adjustment for known scale etc issues.")
        return
    comparison = COMPARISON[key]
    comparison()

if __name__ == "__main__":
    main()
