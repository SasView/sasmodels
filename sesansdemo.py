# Example of conversion of scattering cross section from SANS in absolute
# units into SESANS using a Hankel transformation
# everything is in units of metres except specified otherwise
# Wim Bouwman (w.g.bouwman@tudelft.nl), June 2013

from __future__ import division

from pylab import *
from scipy.special import jv as besselj

# q-range parameters

q = arange(0.0003, 1.0, 0.0003);    # [nm^-1] range wide enough for  Hankel transform
dq=(q[1]-q[0])*1e9;   # [m^-1] step size in q, needed for integration
nq=len(q);
Lambda=2e-10;   # [m] wavelength
# sample parameters
phi=0.1;   # volume fraction
R=100;       # [nm] radius particles
DeltaRho=6e14;  # [m^-2]
V=4/3*pi*R**3 * 1e-27; # [m^3]
th=0.002;    # [m] thickness sample

#2 PHASE SYSTEM
st= 1.5*Lambda**2*DeltaRho**2*th*phi*(1-phi)*R*1e-9  # scattering power in sesans formalism

# Form factor solid sphere
qr=q*R;
P=(3.*(sin(qr)-qr*cos(qr)) / qr**3)**2;
# Structure factor dilute
S=1.;
#2 PHASE SYSTEM
# scattered intensity [m^-1] in absolute units according to SANS 
I=phi*(1-phi)*V*(DeltaRho**2)*P*S;

clf()
subplot(211)  # plot the SANS calculation
plot(q,I,'k')
loglog(q,I)
xlim([0.01, 1])
ylim([1, 1e9])
xlabel(r'$Q [nm^{-1}]$')
ylabel(r'$d\Sigma/d\Omega [m^{-1}]$')


# Hankel transform to nice range for plot
nz=61;
zz=linspace(0,240,nz); # [nm], should be less than reciprocal from q
G=zeros(nz);
for i in range(len(zz)):
    integr=besselj(0,q*zz[i])*I*q;
    G[i]=sum(integr);
G=G*dq*1e9*2*pi; # integr step, conver q into [m**-1] and 2 pi circle integr
# plot(zz,G);
stt= th*Lambda**2/4/pi/pi*G[0]  # scattering power according to SANS formalism
PP=exp(th*Lambda**2/4/pi/pi*(G-G[0]));

subplot(212)
plot(zz,PP,'k',label="Hankel transform") # Hankel transform 1D
xlabel('spin-echo length [nm]')
ylabel('polarisation normalised')
hold(True)

# Cosine transformation of 2D scattering patern
if False:
    qy,qz = meshgrid(q,q)
    qr=R*sqrt(qy**2 + qz**2); # reuse variable names Hankel transform, but now 2D
    P=(3.*(sin(qr)-qr*cos(qr)) / qr**3)**2;
    # Structure factor dilute
    S=1.;
    # scattered intensity [m^-1] in absolute units according to SANS
    I=phi*V*(DeltaRho**2)*P*S;
    GG=zeros(nz);
    for i in range(len(zz)):
        integr=cos(qz*zz[i])*I;
        GG[i]=sum(sum(integr));
    GG=4*GG* dq**2; # take integration step into account take 4 quadrants
    # plot(zz,GG);
    sstt= th*Lambda**2/4/pi/pi*GG[0]  # scattering power according to SANS formalism
    PPP=exp(th*Lambda**2/4/pi/pi*(GG-GG[0]));

    plot(zz,PPP,label="cosine transform") # cosine transform 2D

# For comparison calculation in SESANS formalism, which overlaps perfectly
def gsphere(z,r):
    """
    Calculate SESANS-correlation function for a solid sphere.

    Wim Bouwman after formulae Timofei Kruglov J.Appl.Cryst. 2003 article
    """
    d = z/r
    g = zeros_like(z)
    g[d==0] = 1.
    low = ((d > 0) & (d < 2))
    dlow = d[low]
    dlow2 = dlow**2
    print dlow.shape, dlow2.shape
    g[low] = sqrt(1-dlow2/4.)*(1+dlow2/8.) + dlow2/2.*(1-dlow2/16.)*log(dlow/(2.+sqrt(4.-dlow2)))
    return g

if True:
    plot(zz,exp(st*(gsphere(zz,R)-1)),'r', label="analytical")
legend()
show()
