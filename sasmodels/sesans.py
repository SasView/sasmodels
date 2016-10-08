"""
Conversion of scattering cross section from SANS (I(q), or rather, ds/dO) in absolute
units (cm-1)into SESANS correlation function G using a Hankel transformation, then converting
the SESANS correlation function into polarisation from the SESANS experiment

Everything is in units of metres except specified otherwise (NOT TRUE!!!)
Everything is in conventional units (nm for spin echo length)

Wim Bouwman (w.g.bouwman@tudelft.nl), June 2013
"""

from __future__ import division

import numpy as np  # type: ignore
from numpy import pi, exp  # type: ignore
from scipy.special import jv as besselj
from sas.sasgui.perspectives.fitting.fitpage import FitPage
#import direct_model.DataMixin as model
        
def make_q(q_max, Rmax):
    r"""
    Return a $q$ vector suitable for SESANS covering from $2\pi/ (10 R_{\max})$
    to $q_max$. This is the integration range of the Hankel transform; bigger range and 
    more points makes a better numerical integration.
    Smaller q_min will increase reliable spin echo length range. 
    Rmax is the "radius" of the largest expected object and can be set elsewhere.
    q_max is determined by the acceptance angle of the SESANS instrument.
    """
    from sas.sascalc.data_util.nxsunit import Converter

    q_min = dq = 0.1 * 2*pi / Rmax
    return np.arange(q_min,
                     Converter(q_max[1])(q_max[0],
                                         units="1/A"),
                     dq)
    
def make_all_q(data):
    """
    Return a $q$ vector suitable for calculating the total scattering cross section for
    calculating the effect of finite acceptance angles on Time of Flight SESANS instruments.
    If no acceptance is given, or unwanted (set "unwanted" flag in paramfile), no all_q vector is needed.
    If the instrument has a rectangular acceptance, 2 all_q vectors are needed.
    If the instrument has a circular acceptance, 1 all_q vector is needed
    
    """
    if not data.has_no_finite_acceptance:
        return []
    elif data.has_yz_acceptance(data):
        # compute qx, qy
        Qx, Qy = np.meshgrid(qx, qy)
        return [Qx, Qy]
    else:
        # else only need q
        # data.has_z_acceptance
        return [q]

def transform(data, q_calc, Iq_calc, qmono, Iq_mono):
    """
    Decides which transform type is to be used, based on the experiment data file contents (header)
    (2016-03-19: currently controlled from parameters script)
    nqmono is the number of q vectors to be used for the detector integration
    """
    nqmono = len(qmono)
    #if nqmono == 0 # if radiobutton hankel is active
    if FitPage.hankel.GetValue():
        result = call_hankel(data, q_calc, Iq_calc)
    elif nqmono == 1:
        q = qmono[0]
        result = call_HankelAccept(data, q_calc, Iq_calc, q, Iq_mono)
    else: #if radiobutton Cosine is active
        Qx, Qy = [qmono[0], qmono[1]]
        Qx = np.reshape(Qx, nqx, nqy)
        Qy = np.reshape(Qy, nqx, nqy)
        Iq_mono = np.reshape(Iq_mono, nqx, nqy)
        qx = Qx[0, :]
        qy = Qy[:, 0]
        result = call_Cosine2D(data, q_calc, Iq_calc, qx, qy, Iq_mono)

    return result

def call_hankel(data, q_calc, Iq_calc):
    return hankel((data.x, data.x_unit),
                  (data.lam, data.lam_unit),
                  (data.sample.thickness,
                   data.sample.thickness_unit),
                  q_calc, Iq_calc)
  
def call_HankelAccept(data, q_calc, Iq_calc, q_mono, Iq_mono):
    return hankel(data.x, data.lam * 1e-9,
                  data.sample.thickness / 10,
                  q_calc, Iq_calc)
                  
def call_Cosine2D(data, q_calc, Iq_calc, qx, qy, Iq_mono):
    return hankel(data.x, data.y, data.lam * 1e-9,
                  data.sample.thickness / 10,
                  q_calc, Iq_calc)
                        
def TotalScatter(model, parameters):  #Work in progress!!
#    Calls a model with existing model parameters already in place, then integrate the product of q and I(q) from 0 to (4*pi/lambda)
    allq = np.linspace(0,4*pi/wavelength,1000)
    allIq = 1
    integral = allq*allIq
    


def Cosine2D(wavelength, magfield, thickness, qy, qz, Iqy, Iqz, modelname): #Work in progress!! Needs to call model still
#==============================================================================
#     2D Cosine Transform if "wavelength" is a vector
#==============================================================================
#allq is the q-space needed to create the total scattering cross-section

    Gprime = np.zeros_like(wavelength, 'd')
    s = np.zeros_like(wavelength, 'd')
    sd = np.zeros_like(wavelength, 'd')
    Gprime = np.zeros_like(wavelength, 'd')
    f = np.zeros_like(wavelength, 'd')
    for i, wavelength_i in enumerate(wavelength):
        z = magfield*wavelength_i
        allq=np.linspace() #for calculating the Q-range of the  scattering power integral
        allIq=np.linspace()  # This is the model applied to the allq q-space. Needs to refference the model somehow
        alldq = (allq[1]-allq[0])*1e10
        sigma[i]=wavelength[i]^2*thickness/2/pi*np.sum(allIq*allq*alldq)
        s[i]=1-exp(-sigma)
        for j, Iqy_j, qy_j in enumerate(qy):
            for k, Iqz_k, qz_k in enumerate(qz):
                Iq = np.sqrt(Iqy_j^2+Iqz_k^2)
                q = np.sqrt(qy_j^2 + qz_k^2)
                Gintegral = Iq*cos(z*Qz_k)
                Gprime[i] += Gintegral
#                sigma = wavelength^2*thickness/2/pi* allq[i]*allIq[i]
#                s[i] += 1-exp(Totalscatter(modelname)*thickness)
#                For now, work with standard 2-phase scatter


                sd[i] += Iq
        f[i] = 1-s[i]+sd[i]
        P[i] = (1-sd[i]/f[i])+1/f[i]*Gprime[i]




def HankelAccept(wavelength, magfield, thickness, q, Iq, theta, modelname):
#==============================================================================
#     HankelTransform with fixed circular acceptance angle (circular aperture) for Time of Flight SESANS
#==============================================================================
#acceptq is the q-space needed to create limited acceptance effect
    SElength= wavelength*magfield
    G = np.zeros_like(SElength, 'd')
    threshold=2*pi*theta/wavelength
    for i, SElength_i in enumerate(SElength):
        allq=np.linspace() #for calculating the Q-range of the  scattering power integral
        allIq=np.linspace()  # This is the model applied to the allq q-space. Needs to refference the model somehow
        alldq = (allq[1]-allq[0])*1e10
        sigma[i]=wavelength[i]^2*thickness/2/pi*np.sum(allIq*allq*alldq)
        s[i]=1-exp(-sigma)

        dq = (q[1]-q[0])*1e10
        a = (x<threshold)
        acceptq = a*q
        acceptIq = a*Iq

        G[i] = np.sum(besselj(0, acceptq*SElength_i)*acceptIq*acceptq*dq)

#        G[i]=np.sum(integral)

    G *= dq*1e10*2*pi

    P = exp(thickness*wavelength**2/(4*pi**2)*(G-G[0]))
    
def hankel(SElength, wavelength, thickness, q, Iq):
    r"""
    Compute the expected SESANS polarization for a given SANS pattern.

    Uses the hankel transform followed by the exponential.  The values for *zz*
    (or spin echo length, or delta), wavelength and sample thickness should
    come from the dataset.  $q$ should be chosen such that the oscillations
    in $I(q)$ are well sampled (e.g., $5 \cdot 2 \pi/d_{\max}$).

    *SElength* [A] is the set of $z$ points at which to compute the
    Hankel transform

    *wavelength* [m]  is the wavelength of each individual point *zz*

    *thickness* [cm] is the sample thickness.

    *q* [A$^{-1}$] is the set of $q$ points at which the model has been
    computed. These should be equally spaced.

    *I* [cm$^{-1}$] is the value of the SANS model at *q*
    """

    from sas.sascalc.data_util.nxsunit import Converter
    wavelength = Converter(wavelength[1])(wavelength[0],"A")
    thickness = Converter(thickness[1])(thickness[0],"A")
    Iq = Converter("1/cm")(Iq,"1/A") # All models default to inverse centimeters
    SElength = Converter(SElength[1])(SElength[0],"A")

    G = np.zeros_like(SElength, 'd')
#==============================================================================
#     Hankel Transform method if "wavelength" is a scalar; mono-chromatic SESANS
#==============================================================================
    for i, SElength_i in enumerate(SElength):
        integral = besselj(0, q*SElength_i)*Iq*q
        G[i] = np.sum(integral)
    G0 = np.sum(Iq*q)

    # [m^-1] step size in q, needed for integration
    dq = (q[1]-q[0])

    # integration step, convert q into [m**-1] and 2 pi circle integration
    G *= dq*2*pi
    G0 = np.sum(Iq*q)*dq*2*np.pi

    P = exp(thickness*wavelength**2/(4*pi**2)*(G-G0))

    return P
