from __future__ import division, print_function
# Make sasmodels available on the path
import sys,os
BETA_DIR = os.path.dirname(os.path.realpath(__file__))
SASMODELS_DIR = os.path.dirname(os.path.dirname(BETA_DIR))
sys.path.insert(0, SASMODELS_DIR)
import os
from matplotlib import pyplot as plt
import numpy as np
from numpy import pi, sin, cos, sqrt, fabs
from numpy.polynomial.legendre import leggauss
from scipy.special import j1 as J1
from numpy import inf
from scipy.special import gammaln  # type: ignore

# THE FOLLOWING 7 BOOLEANS TAILOR WHICH GRAPHS ARE PRINTED WHEN THE PROGRAM RUNS
#RICHARD, YOU WILL WANT SASVIEW, SASFIT, SCHULZ, ELLIPSOID, AND SPHERE ALL TRUE.
ELLIPSOID = False
SPHERE = True

GAUSSIAN = True
SCHULZ = False

SASVIEW=True
SASFIT=True
YUN = True

def data_file(name):
    return os.path.join(BETA_DIR + '\\data_files', name)

#Used to calculate F(q) for the cylinder, sphere, ellipsoid models
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
        c2 =  1.0/280.0
        c3 = -1.0/15120.0
        c4 =  1.0/1330560.0
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
    calculator.pars.setdefault('background', 1e-3)
    return calculator

#gives the hardsphere structure factor that sasview uses
def hardsphere_simple(q, radius_effective, volfraction): 
    CUTOFFHS=0.05 
    if fabs(radius_effective) < 1.E-12:
        HARDSPH=1.0
        return HARDSPH
    X = 1.0/( 1.0 -volfraction)
    D= X*X
    A= (1.+2.*volfraction)*D
    A *=A
    X=fabs(q*radius_effective*2.0)
    if X < 5.E-06:
        HARDSPH=1./A
        return HARDSPH
    X2 =X*X
    B = (1.0 +0.5*volfraction)*D
    B *= B
    B *= -6.*volfraction
    G=0.5*volfraction*A
    if X < CUTOFFHS:
        FF = 8.0*A +6.0*B + 4.0*G + ( -0.8*A -B/1.5 -0.5*G +(A/35. +0.0125*B +0.02*G)*X2)*X2
        HARDSPH= 1./(1. + volfraction*FF )
        return HARDSPH   
    X4=X2*X2
    S, C = sin(X), cos(X)
    FF=  (( G*( (4.*X2 -24.)*X*S -(X4 -12.*X2 +24.)*C +24. )/X2 + B*(2.*X*S -(X2-2.)*C -2.) )/X + A*(S-X*C))/X ;
    HARDSPH= 1./(1. + 24.*volfraction*FF/X2 );
    return HARDSPH

#Used in gaussian quadrature for polydispersity
#returns values and the probability of those values based on gaussian distribution
def gaussian_distribution(center, sigma,lb,ub):
    #3 standard deviations covers approx. 99.7% 
    if sigma != 0:
        nsigmas=3
        x = np.linspace(center-sigma*nsigmas, center+sigma*nsigmas, num=35)
        x= x[(x >= lb) & (x <= ub)]
        px = np.exp((x-center)**2 / (-2.0 * sigma * sigma))
        return x, px
    else:
        return np.array([center]), np.array([1])

def schulz_distribution(center, sigma, lb, ub):
    if sigma != 0:
        nsigmas=8
        x = np.linspace(center-sigma*nsigmas, center+sigma*nsigmas, num=80)
        x= x[(x >= lb) & (x <= ub)]
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
    if (radius_polar * radius_equatorial != 0):
        bd = 1.0 - ee
        e1 = np.sqrt(ee)
        b1 = 1.0 + np.arcsin(e1) / (e1*np.sqrt(bd))
        bL = (1.0 + e1) / (1.0 - e1)
        b2 = 1.0 + bd / 2 / e1 * np.log(bL)
        delta = 0.75 * b1 * b2
    ddd = np.zeros_like(radius_polar)
    ddd = 2.0*(delta + 1.0)*radius_polar*radius_equatorial**2
    return 0.5*ddd**(1.0 / 3.0)

def ellipsoid_volume(radius_polar,radius_equatorial):
    volume = (4./3.)*pi*radius_polar*radius_equatorial**2
    return volume

# F1 is F(q)
# F2 is F(g)^2
#IQM is I(q) with monodispersity
#IQSM is I(q) with structure factor S(q) and monodispersity
#IQBM is I(q) with Beta Approximation and monodispersity
#SQ is monodisperse approach for structure factor
#SQ_EFF is the effective structure factor from beta approx
def ellipsoid_Theta(q, radius_polar, radius_equatorial, sld, sld_solvent,volfraction=0,effective_radius=0):
    #creates values z and corresponding probabilities w from legendre-gauss quadrature
    z, w = leggauss(76)
    F1 = np.zeros_like(q)
    F2 = np.zeros_like(q)
    #use a u subsition(u=cos) and then u=(z+1)/2 to change integration from
    #0->2pi with respect to alpha to -1->1 with respect to z, allowing us to use 
    #legendre-gauss quadrature
    for k, qk in enumerate(q):
        r = sqrt(radius_equatorial**2*(1-((z+1)/2)**2)+radius_polar**2*((z+1)/2)**2)
        F2i = ((sld-sld_solvent)*ellipsoid_volume(radius_polar,radius_equatorial)*sas_3j1x_x(qk*r))**2
        F2[k] = np.sum(w*F2i)
        F1i = (sld-sld_solvent)*ellipsoid_volume(radius_polar,radius_equatorial)*sas_3j1x_x(qk*r)
        F1[k] = np.sum(w*F1i)
    #the 1/2 comes from the change of variables mentioned above
    F2 = F2/2.0
    F1 = F1/2.0
    if effective_radius == 0:
        effective_radius = ER_ellipsoid(radius_polar,radius_equatorial)
    else:
        effective_radius = effective_radius
    SQ = np.array([hardsphere_simple(qk, effective_radius, volfraction) for qk in q])   
    SQ_EFF = 1 + F1**2/F2*(SQ - 1) 
    IQM = 1e-4/ellipsoid_volume(radius_polar,radius_equatorial)*F2
    IQSM = volfraction*IQM*SQ
    IQBM = volfraction*IQM*SQ_EFF
    return F1, F2, IQM, IQSM, IQBM, SQ, SQ_EFF

#IQD is I(q) polydispursed, IQSD is I(q)S(q) polydispursed, etc. 
#IQBD HAS NOT BEEN CROSS CHECKED AT ALL
def ellipsoid_pe(q, Rp, Re, sld, sld_solvent, volfraction=0,effective_radius=0,radius_polar_pd=0.1,radius_equatorial_pd=0.1,distribution='gaussian'):
    if distribution == 'gaussian':
        Rp_val, Rp_prob = gaussian_distribution(Rp, radius_polar_pd*Rp, 0, inf)
        Re_val, Re_prob = gaussian_distribution(Re, radius_equatorial_pd*Re, 0, inf)
    elif distribution == 'schulz':
        Rp_val, Rp_prob = schulz_distribution(Rp, radius_polar_pd*Rp, 0, inf)
        Re_val, Re_prob = schulz_distribution(Re, radius_equatorial_pd*Re, 0, inf)
    Normalization = 0
    total_weight = 0
    PQ = np.zeros_like(q)
    F12,F21 = np.zeros_like(q), np.zeros_like(q)
    radius_eff = 0
    for k, Rpk in enumerate(Rp_val):
        for i, Rei in enumerate(Re_val):
            F1i, F2i, PQi, IQSM, IQBM, SQ, SQ_EFF = ellipsoid_Theta(q,Rpk,Rei,sld,sld_solvent)
            radius_eff += Rp_prob[k]*Re_prob[i]*ER_ellipsoid(Rpk,Rei)
            total_weight +=  Rp_prob[k]*Re_prob[i]
            Normalization += Rp_prob[k]*Re_prob[i]*ellipsoid_volume(Rpk, Rei)
            PQ += PQi*Rp_prob[k]*Re_prob[i]*ellipsoid_volume(Rpk, Rei)
            F12 += F1i*Rp_prob[k]*Re_prob[i]*ellipsoid_volume(Rpk, Rei)
            F21 += F2i*Rp_prob[k]*Re_prob[i]*ellipsoid_volume(Rpk, Rei)
    F12 = (F12/Normalization)**2
    F21 = F21/Normalization
    if effective_radius == 0:
        effective_radius = radius_eff/total_weight
    else:
        effective_radius = effective_radius
    SQ = np.array([hardsphere_simple(qk, effective_radius, volfraction) for qk in q])   
    SQ_EFF = 1 + F12/F21*(SQ - 1)     
    IQD = PQ/Normalization
    IQSD = volfraction*IQD*SQ
    IQBD = volfraction*IQD*SQ_EFF
    return IQD, IQSD, IQBD, SQ, SQ_EFF

#polydispersity for sphere
def sphere_r(q,radius,sld,sld_solvent,volfraction=0,radius_pd=0.1,distribution='gaussian',norm_type='volfraction'):
    if distribution == 'gaussian':    
        radius_val, radius_prob = gaussian_distribution(radius, radius_pd*radius, 0, inf)
    elif distribution == 'schulz':
        radius_val, radius_prob = schulz_distribution(radius, radius_pd*radius, 0, inf)
    Normalization=0
    F2 = np.zeros_like(q)
    F1 = np.zeros_like(q)
    if norm_type == 'numdensity':
        for k, rk in enumerate(radius_val):
            volume = 4./3.*pi*rk**3
            Normalization += radius_prob[k]*volume
            F2i = ((sld-sld_solvent)*volume*sas_3j1x_x(q*rk))**2/volume
            F1i = (sld-sld_solvent)*volume*sas_3j1x_x(q*rk)/volume
            
            F2 += radius_prob[k]*volume*F2i
            F1 += radius_prob[k]*volume*F1i
        F21 = F2/Normalization
        F12 = (F1/Normalization)**2
        SQ = np.array([hardsphere_simple(qk, radius, volfraction) for qk in q])
        SQ_EFF = 1 + F12/F21*(SQ - 1) 
        IQD = 1e-4*F21
        IQSD = volfraction*IQD*SQ
        IQBD = volfraction*IQD*SQ_EFF
    elif norm_type == 'volfraction':
        for k, rk in enumerate(radius_val):
            volume = 4./3.*pi*rk**3
            Normalization += radius_prob[k]
            F2i = ((sld-sld_solvent)*volume*sas_3j1x_x(q*rk))**2
            F1i = (sld-sld_solvent)*volume*sas_3j1x_x(q*rk)
            F2 += radius_prob[k]*F2i
            F1 += radius_prob[k]*F1i
    
        F21 = F2/Normalization
        F12 = (F1/Normalization)**2
        SQ = np.array([hardsphere_simple(qk, radius, volfraction) for qk in q])
        SQ_EFF = 1 + F12/F21*(SQ - 1) 
        IQD = 1e-4/(4./3.*pi*radius**3)*F21
        IQSD = volfraction*IQD*SQ
        IQBD = volfraction*IQD*SQ_EFF
    return IQD, IQSD, IQBD, SQ, SQ_EFF

###############################################################################
###############################################################################
###############################################################################
##################                                           ##################
##################                   TESTS                   ##################
##################                                           ##################
###############################################################################
###############################################################################
###############################################################################
# SASVIEW    
if (ELLIPSOID == True) & (GAUSSIAN == True) & (SASVIEW == True):
    q = np.logspace(-5, 0, 50)
    volfraction = 0.2
    model = build_model("ellipsoid",q)
    Sq = model(radius_polar=20,radius_equatorial=400,sld=4,sld_solvent=1,\
               background=0,radius_polar_pd=.1,radius_polar_pd_n=35,radius_equatorial_pd=.1,radius_equatorial_pd_n=35)
    IQD, IQSD, IQBD, SQ, SQ_EFF = ellipsoid_pe(q,20,400,4,1,volfraction)
    Sq2 = model(radius_polar=20,radius_equatorial=10,sld=2,sld_solvent=1,background=0)
    IQM = ellipsoid_Theta(q,20,10,2,1)[2]
    model2 = build_model("ellipsoid@hardsphere", q)
    Sq3 = model2(radius_polar=20,radius_equatorial=400,sld=4,sld_solvent=1,background=0,radius_polar_pd=.1,\
                 radius_polar_pd_n=35,radius_equatorial_pd=.1,radius_equatorial_pd_n=35)
    Sqp = model(radius_polar=20,radius_equatorial=10,sld=4,sld_solvent=1,\
                   background=0,radius_polar_pd=0.1,radius_polar_pd_n=35,radius_equatorial_pd=0,radius_equatorial_pd_n=35)
    IQDp = ellipsoid_pe(q,20,10,4,1,radius_polar_pd=0.1,radius_equatorial_pd=0)[0]  
    Sqp2 = model(radius_polar=20,radius_equatorial=10,sld=4,sld_solvent=1,\
                   background=0,radius_polar_pd=0,radius_polar_pd_n=35,radius_equatorial_pd=0.1,radius_equatorial_pd_n=35)
    IQDe = ellipsoid_pe(q,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0.1)[0]
    print('\n\n\n\n ELLIPSOID COMPARISON WITH SASVIEW WITHOUT V')
    plt.subplot(323)
    plt.loglog(q, Sq,'r')
    plt.loglog(q,IQD,'b')
    plt.title('IQD')
    plt.subplot(324)
    plt.semilogx(q,Sq/IQD-1)
    plt.subplot(321)
    plt.title('IQM')
    plt.loglog(q, Sq2,'r')
    plt.loglog(q,IQM,'b')
    plt.subplot(322)
    plt.semilogx(q,IQM/Sq2-1)
    plt.subplot(325)
    plt.title('IQSD')
    plt.loglog(q, Sq3, '-b')
    plt.loglog(q, IQSD, '-r')
    plt.subplot(326)
    plt.semilogx(q, Sq3/IQSD-1, 'b')
    plt.tight_layout()
    plt.show()
    plt.subplot(221)
    plt.loglog(q,IQDp,'r')
    plt.loglog(q,Sqp,'b')
    plt.title('IQD Polar')
    plt.subplot(222)
    plt.plot(q,IQDp/Sqp-1)
    plt.subplot(223)
    plt.loglog(q,IQDe,'r')
    plt.loglog(q,Sqp2,'b')
    plt.title('IQD Equatorial')
    plt.subplot(224)
    plt.plot(q,IQDe/Sqp2-1)
    plt.tight_layout()
    plt.show()
#comparing monodisperse Beta approximation to version without
    q=np.linspace(1e-5,1)
    F1, F2, IQM, IQSM, IQBM, SQ, SQ_EFF = ellipsoid_Theta(q,20,10,2,1,0.15,12.59921049894873)
    plt.subplot(321)
    plt.loglog(q,IQBM,'g',label='IQBM')
    plt.loglog(q,IQSM,'r',label= 'IQSM')
    plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.subplot(323)
    plt.plot(q,IQBM,'g',label='IQBM')
    plt.plot(q,IQSM,'r',label= 'IQSM')
    plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.subplot(325)
    plt.plot(q,IQBM/IQSM,'r',label= 'IQBM/IQSM')
    plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.tight_layout()
    plt.show()
#comparing polydispersed Beta approximation to version without
    IQD, IQSD, IQBD, SQ, SQ_EFF = ellipsoid_pe(q, 20,10, 2,1, 0.15,12.59921049894873)
    plt.subplot(321)
    plt.loglog(q, SQ)
    plt.loglog(q, SQ_EFF,label='SQ, SQ_EFF')
    plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.subplot(323)
    plt.loglog(q,IQBD)
    plt.loglog(q,IQSD,label='IQBD,IQSD')
    plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.subplot(325)
    plt.plot(q,IQBD)
    plt.plot(q,IQSD,label='IQBD,IQSD')
    plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.tight_layout()
    plt.show()

# YUN
if (ELLIPSOID == True) & (GAUSSIAN == True) & (YUN == True):
    #Yun's data for Beta Approximation
    data = np.loadtxt(data_file('yun_ellipsoid.dat'),skiprows=2).T   
    print('\n\n ELLIPSOID COMPARISON WITH YUN WITHOUT V')
    Q=data[0]
    F1,F2,IQM,IQSM,IQBM,SQ,SQ_EFF = ellipsoid_Theta(Q,20,10,2,1,0.15,12.59921049894873)
    plt.subplot(211)
    plt.loglog(Q,0.15*data[3]*data[5]*ellipsoid_volume(20,10)*1e-4)
    plt.loglog(Q,IQSM,'-g',label='IQSM')
    plt.loglog(data[0], data[3]*ellipsoid_volume(20,10)*1e-4,'-k')
    plt.loglog(Q,IQM,'r',label = 'IQM')
    plt.ylim(1e-5,4)
    plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.show()
    plt.subplot(221)
    plt.loglog(data[0],SQ)
    plt.loglog(data[0],data[5])
    plt.title('SQ')
    plt.subplot(222)
    plt.semilogx(data[0],SQ/data[5]-1)
    plt.subplot(223)
    plt.title('SQ_EFF')
    plt.loglog(data[0],SQ_EFF)
    plt.loglog(data[0],data[6])
    plt.subplot(224)
    plt.semilogx(data[0],(SQ_EFF/data[6])-1,label='Sq effective ratio')
    plt.tight_layout()
    plt.show()

#SASFIT
if ELLIPSOID and GAUSSIAN and SASFIT:
    print('\n\n ELLIPSOID COMPARISON WITH SASFIT WITHOUT V')
    #Rp=20,Re=10,eta_core=4,eta_solv=1
    sasdata = np.loadtxt(data_file('sasfit_ellipsoid_IQM.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(321)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0)[0],'b')
    plt.loglog(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata),'g')
    plt.title('IQM')
    plt.subplot(322)
    plt.plot(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0)[0]-1)
    plt.tight_layout()
    plt.show()
    print('========================================================================================')
    print('========================================================================================')
#N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.15
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(321)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.1,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 10% 0.15')
    plt.subplot(322)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.1,radius_equatorial_pd=0)[3]-1)
#N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.15
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq2.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(323)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.3,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 30% .15')
    plt.subplot(324)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.3,radius_equatorial_pd=0)[3]-1)
#N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.15
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq3.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(325)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.6,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 60% .15')
    plt.subplot(326)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.6,radius_equatorial_pd=0)[3]-1)
    plt.tight_layout()
    plt.show() 
    print('========================================================================================')
    print('========================================================================================')
#N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.3
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq4.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(321)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.1,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 10% .3')
    plt.subplot(322)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.1,radius_equatorial_pd=0)[3]-1)
#N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.3
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq5.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(323)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.3,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 30% .3')
    plt.subplot(324)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.3,radius_equatorial_pd=0)[3]-1)
#N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.3
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq6.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(325)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.6,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 60% .3')
    plt.subplot(326)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.6,radius_equatorial_pd=0)[3]-1)
    plt.tight_layout()
    plt.show()
    print('========================================================================================')
    print('========================================================================================')
#N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.6
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq7.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(321)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.1,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 10% .6')
    plt.subplot(322)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.1,radius_equatorial_pd=0)[3]-1)
#N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.6
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq8.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(323)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.3,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 30% .6')
    plt.subplot(324)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.3,radius_equatorial_pd=0)[3]-1)
#N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.6
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sq9.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(325)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.6,radius_equatorial_pd=0)[3])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ poly 60% .6')
    plt.subplot(326)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.6,radius_equatorial_pd=0)[3]-1)
    plt.tight_layout()
    plt.show()
    print('========================================================================================')
    print('========================================================================================')
#checks beta structre factor
#N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.15
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(321)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.1,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 10% 0.15')
    plt.subplot(322)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.1,radius_equatorial_pd=0)[4]-1)
#N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.15
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff2.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(323)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.3,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 30% .15')
    plt.subplot(324)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.3,radius_equatorial_pd=0)[4]-1)
#N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.15
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff3.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(325)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.6,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 60% .15')
    plt.subplot(326)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,.15,radius_polar_pd=0.6,radius_equatorial_pd=0)[4]-1)
    plt.tight_layout()
    plt.show() 
    print('========================================================================================')
    print('========================================================================================')
#N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.3
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff4.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(321)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.1,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 10% .3')
    plt.subplot(322)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.1,radius_equatorial_pd=0)[4]-1)
#N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.3
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff5.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(323)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.3,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 30% .3')
    plt.subplot(324)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.3,radius_equatorial_pd=0)[4]-1)
#N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.3
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff6.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(325)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.6,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 60% .3')
    plt.subplot(326)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.3,radius_polar_pd=0.6,radius_equatorial_pd=0)[4]-1)
    plt.tight_layout()
    plt.show()
    print('========================================================================================')
    print('========================================================================================')
#N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.1354236254,.6
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff7.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(321)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.1,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 10% .6')
    plt.subplot(322)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.1,radius_equatorial_pd=0)[4]-1)
#N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.0901197149,.6
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff8.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(323)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.3,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 30% .6')
    plt.subplot(324)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.3,radius_equatorial_pd=0)[4]-1)
#N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, hardsphere ,13.336060917,.6
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_ellipsoid_sqeff9.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(325)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.6,radius_equatorial_pd=0)[4])
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.title('SQ_EFF poly 60% .6')
    plt.subplot(326)
    plt.plot(sasqdata,np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,0.6,radius_polar_pd=0.6,radius_equatorial_pd=0)[4]-1)
    plt.tight_layout()
    plt.show()
    print('========================================================================================')
    print('========================================================================================')
#checks polydispersion for single parameter and varying pd with sasfit
#N=1,s=2,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    sasdata = np.loadtxt(data_file('sasfit_ellipsoid_IQD.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(221)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0.1,radius_equatorial_pd=0)[0])
    plt.loglog(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata))
    plt.title('IQD polar 10%')
    plt.subplot(222)
    plt.plot(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0.1,radius_equatorial_pd=0)[0]-1)
#N=1,s=1,X0=10,distr equatorial Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    sasdata = np.loadtxt(data_file('sasfit_ellipsoid_IQD2.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(223)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0.1)[0])
    plt.loglog(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata))
    plt.title('IQD equatorial 10%')
    plt.subplot(224)
    plt.plot(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0.1)[0]-1)
    plt.tight_layout()
    plt.show()
    print('========================================================================================')
    print('========================================================================================')
    #N=1,s=6,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    sasdata = np.loadtxt(data_file('sasfit_ellipsoid_IQD3.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(221)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0.3,radius_equatorial_pd=0)[0])
    plt.loglog(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata))
    plt.title('IQD polar 30%')
    plt.subplot(222)
    plt.plot(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0.3,radius_equatorial_pd=0)[0]-1)
#N=1,s=3,X0=10,distr equatorial Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    sasdata = np.loadtxt(data_file('sasfit_ellipsoid_IQD4.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(223)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0.3)[0])
    plt.loglog(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata))
    plt.title('IQD equatorial 30%')
    plt.subplot(224)
    plt.plot(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0.3)[0]-1)
    plt.tight_layout()
    plt.show()
    print('========================================================================================')
    print('========================================================================================')
    #N=1,s=12,X0=20,distr polar Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    sasdata = np.loadtxt(data_file('sasfit_ellipsoid_IQD5.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(221)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0.6,radius_equatorial_pd=0)[0])
    plt.loglog(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata))
    plt.title('IQD polar 60%')
    plt.subplot(222)
    plt.plot(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0.6,radius_equatorial_pd=0)[0]-1)
#N=1,s=6,X0=10,distr equatorial Rp=20,Re=10,eta_core=4,eta_solv=1, no structure poly
    sasdata = np.loadtxt(data_file('sasfit_ellipsoid_IQD6.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(223)
    plt.loglog(sasqdata, ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0.6)[0])
    plt.loglog(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata))
    plt.title('IQD equatorial 60%')
    plt.subplot(224)
    plt.plot(sasqdata,1e-4/ellipsoid_volume(20,10)*np.array(sasvaluedata)/ellipsoid_pe(sasqdata,20,10,4,1,radius_polar_pd=0,radius_equatorial_pd=0.6)[0]-1)
    plt.tight_layout()
    plt.show()
    print('========================================================================================')
    print('========================================================================================')

# TEST FOR RICHARD
    #radius=20,sld=4,sld_solvent=1,volfraction=0.3,radius_pd=0.1
    #We have scaled the output from sasfit by 1e-4*volume*volfraction
    #0.10050378152592121
if (SPHERE == True) & (SCHULZ == True) & (SASVIEW == True) & (SASFIT == True):
    print('                   *****SPHERE MODEL*****')
    sasdata = np.loadtxt(data_file('richard_test.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    model=build_model('sphere',sasqdata)
    Iq = model(radius=20,radius_pd=0.1,radius_pd_n=80,sld=4,sld_solvent=1,background=0,radius_pd_type='schulz',radius_pd_nsigma=8)
    model2=build_model('sphere@hardsphere',sasqdata)
    IQS = model2(radius=20,radius_pd=0.1,radius_pd_n=80,sld=4,sld_solvent=1,background=0,radius_pd_type='schulz',radius_pd_nsigma=8,volfraction=0.3)
    
    sasvaluedata=map(float,sasdata[1])
    sasdata2 = np.loadtxt(data_file('richard_test2.txt'),dtype=str,delimiter=';').T
    sasqdata2=map(float,sasdata2[0])
    sasvaluedata2=map(float,sasdata2[1])
    sasdata3 = np.loadtxt(data_file('richard_test3.txt'),dtype=str,delimiter=';').T
    sasqdata3=map(float,sasdata3[0])
    sasvaluedata3=map(float,sasdata3[1])
    IQD, IQSD, IQBD, SQ, SQ_EFF = sphere_r(np.array(sasqdata),20,4,1,.3,0.1,distribution='schulz')
    plt.grid(True)
    plt.title('IQD(blue) vs. SASVIEW(red) vs. SASFIT(green)')
    plt.loglog(sasqdata,Iq,'b')
    plt.loglog(sasqdata,IQD,'r')
    plt.loglog(sasqdata,1e-4/(4./3.*pi*20**3)*np.array(sasvaluedata),'g')
    plt.show()
    plt.grid(True)
    plt.title('IQSD(blue) vs. SASVIEW(red) vs. SASFIT(green)')
    plt.loglog(sasqdata,IQSD,'b')
    plt.loglog(sasqdata,IQS,'r')
    plt.loglog(sasqdata,1e-4/(4./3.*pi*20**3)*0.3*np.array(sasvaluedata2),'g')
    plt.show()
    plt.grid(True)
    plt.title('IQBD(blue) vs. SASFIT(green)')
    plt.loglog(sasqdata,IQBD,'b')
    plt.loglog(sasqdata,1e-4/(4./3.*pi*20**3)*0.3*np.array(sasvaluedata3),'g')
    plt.show()

# TEST FOR RICHARD
if (ELLIPSOID == True) & (SCHULZ == True) & (SASVIEW == True) & (SASFIT == True):   
    #polarradius=20, equatorialradius=10, sld=4,sld_solvent=1,volfraction=0.3,radius_polar_pd=0.1
         #Effective radius =13.1353356684
        #We have scaled the output from sasfit by 1e-4*volume*volfraction
        #0.10050378152592121
    print('                   *****ELLIPSOID MODEL*****')
    sasdata = np.loadtxt(data_file('richard_test4.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    model=build_model('ellipsoid',sasqdata)
    Iq = model(radius_polar=20,radius_polar_pd=0.1,radius_polar_pd_n=80, radius_equatorial=10, radius_equatorial_pd=0, sld=4,sld_solvent=1,background=0,radius_polar_pd_type='schulz',radius_polar_pd_nsigma=8)
    model2=build_model('ellipsoid@hardsphere',sasqdata)
    IQS =  model2(radius_polar=20,radius_equatorial=10,sld=4,sld_solvent=1,background=0,radius_polar_pd=.1,\
                     radius_polar_pd_n=80,radius_equatorial_pd=0,radius_polar_pd_type='schulz',radius_polar_pd_nsigma=8,volfraction=0.3)
    sasvaluedata=map(float,sasdata[1])
    sasdata2 = np.loadtxt(data_file('richard_test5.txt'),dtype=str,delimiter=';').T
    sasqdata2=map(float,sasdata2[0])
    sasvaluedata2=map(float,sasdata2[1])
    sasdata3 = np.loadtxt(data_file('richard_test6.txt'),dtype=str,delimiter=';').T
    sasqdata3=map(float,sasdata3[0])
    sasvaluedata3=map(float,sasdata3[1])
    IQD, IQSD, IQBD, SQ, SQ_EFF = ellipsoid_pe(np.array(sasqdata),20,10,4,1,.3,radius_polar_pd=0.1,radius_equatorial_pd=0,distribution='schulz')    
    plt.grid(True)
    plt.title('IQD(blue) vs. SASVIEW(red) vs. SASFIT(green)')
    plt.loglog(sasqdata,Iq,'b')
    plt.loglog(sasqdata,IQD,'r')
    plt.loglog(sasqdata,1e-4/(4./3.*pi*20*10**2)*np.array(sasvaluedata),'g')
    plt.show()
    plt.grid(True)
    plt.title('IQSD(blue) vs. SASVIEW(red) vs. SASFIT(green)')
    plt.loglog(sasqdata,IQSD,'b')
    plt.loglog(sasqdata,IQS,'r')
    plt.loglog(sasqdata,1e-4/(4./3.*pi*20*10**2)*0.3*np.array(sasvaluedata2),'g')
    plt.show()
    plt.grid(True)
    plt.title('IQBD(blue) vs. SASFIT(green)')
    plt.loglog(sasqdata,IQBD,'b')
    plt.loglog(sasqdata,1e-4/(4./3.*pi*20*10**2)*0.3*np.array(sasvaluedata3),'g')
    plt.show()    
    
# SASVIEW GAUSS
if (SPHERE == True) & (GAUSSIAN == True) & (SASVIEW == True):
    q = np.logspace(-5, 0, 200)
    IQD, IQSD, IQBD, SQ, SQ_EFF = sphere_r(q,20,4,1,.3)
    model = build_model("sphere", q)
    Sq = model(radius=20,radius_pd=0.1,radius_pd_n=35,sld=4,sld_solvent=1,background=0)
    model2 = build_model("sphere@hardsphere", q)
    S=build_model("hardsphere",q)(radius_effective=20,volfraction=.3,background=0)
    Sq2 = model2(radius=20,radius_pd=0.1,radius_pd_n=35,sld=4,sld_solvent=1,background=0,volfraction=.3)
    print('\n\n SPHERE COMPARISON WITH SASVIEW WITHOUT V')
    plt.subplot(221)
    plt.title('IQD')
    plt.loglog(q, IQD, '-b')
    plt.loglog(q, Sq, '-r')
    plt.subplot(222)
    plt.semilogx(q, Sq/IQD-1, '-g')
    plt.tight_layout()
    plt.show()
    plt.subplot(221)
    plt.title('SQ')
    plt.plot(q, SQ, '-r')
    plt.plot(q,S,'-k')
    plt.subplot(222)
    plt.plot(SQ/S-1)
    plt.tight_layout()
    plt.show()
    plt.subplot(221)
    plt.title('IQSD')
    plt.loglog(q, IQSD, '-b')
    plt.loglog(q, Sq2, '-r',label='IQSD')
    plt.subplot(222)
    plt.semilogx(q, Sq2/IQSD-1, '-g',label='IQSD ratio')
    plt.tight_layout()
    plt.show()
    IQD, IQSD, IQBD, SQ, SQ_EFF = sphere_r(q,20,4,1,0.15)
    plt.subplot(211)
    plt.title('SQ vs SQ_EFF')
    plt.plot(q, SQ)
    plt.plot(q, SQ_EFF)
    plt.tight_layout()
    plt.show() 
    plt.subplot(221)
    plt.title('IQSD vs IQBD')
    plt.loglog(q,IQBD)
    plt.loglog(q,IQSD,label='IQBD,IQSD')
    plt.subplot(222)
    plt.plot(q,IQBD)
    plt.plot(q,IQSD,)
    plt.tight_layout()
    plt.show() 

# SASFIT GAUSS
if (SPHERE == True) & (GAUSSIAN == True) & (SASFIT == True):
    print('\n\n SPHERE COMPARISON WITH SASFIT WITHOUT V')
#N=1,s=2,X0=20,distr radius R=20,eta_core=4,eta_solv=1,.3
    sasdata = np.loadtxt(data_file('sasfit_sphere_IQD.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    IQD, IQSD, IQBD, SQ, SQ_EFF = sphere_r(np.array(sasqdata),20,4,1,.3)
    plt.subplot(221)
    plt.title('IQD')
    plt.loglog(sasqdata,IQD )
    plt.loglog(sasqdata,1e-4/(4./3*pi*20**3)*np.array(sasvaluedata))
    plt.subplot(222)
    plt.semilogx(sasqdata,1e-4/(4./3*pi*20**3)*np.array(sasvaluedata)/IQD-1)
    plt.tight_layout()
    plt.show()
    sasdata = np.loadtxt(data_file('sasfit_sphere_IQSD.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(221)
    plt.title('IQSD')
    plt.loglog(sasqdata, IQSD)
    plt.loglog(sasqdata,1e-4/(4./3*pi*20**3)*0.3*np.array(sasvaluedata))
    plt.subplot(222)
    plt.semilogx(sasqdata,1e-4/(4./3*pi*20**3)*0.3*np.array(sasvaluedata)/IQSD-1)
    plt.tight_layout()
    plt.show()
    sasdata = np.loadtxt(data_file('sasfit_sphere_IQBD.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(221)
    plt.title('IQBD')
    plt.loglog(sasqdata,IQBD)
    plt.loglog(sasqdata,1e-4/(4./3*pi*20**3)*0.3*np.array(sasvaluedata))
    plt.subplot(222)
    plt.semilogx(sasqdata,1e-4/(4./3*pi*20**3)*0.3*np.array(sasvaluedata)/IQBD-1)
    plt.tight_layout()
    plt.show()
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_sphere_sq.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(221)
    plt.title('SQ')
    plt.loglog(sasqdata, SQ)
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.subplot(222)
    plt.semilogx(sasqdata,np.array(sasvaluedata)/SQ-1)
    plt.tight_layout()
    plt.show()
    sasdata = np.loadtxt(data_file('sasfit_polydisperse_sphere_sqeff.txt'),dtype=str,delimiter=';').T
    sasqdata=map(float,sasdata[0])
    sasvaluedata=map(float,sasdata[1])
    plt.subplot(221)
    plt.title('SQ_EFF')
    plt.loglog(sasqdata,SQ_EFF)
    plt.loglog(sasqdata,np.array(sasvaluedata))
    plt.subplot(222)
    plt.semilogx(sasqdata,np.array(sasvaluedata)/SQ_EFF-1)
    plt.tight_layout()
    plt.show()