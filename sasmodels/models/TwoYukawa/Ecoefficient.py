from typing import Callable, Tuple

import numpy as np
from numpy import exp, pi, cos, sin, cosh
from numpy.typing import NDArray

# CalCoeff.m
# CheckList.m
# Ck, Cr, Cx

class TYCoeff:
    phi: float
    z: Tuple[float, float]
    bigK: Tuple[float, float]
    k: Tuple[float, float]
    # Unused
    d1Factor: float = 1.
    d2Factor: float = 1.

    Ecoefficient: NDArray
    ABC1C2: # Callable[[float, float], [float, float, float, float]]

    def __init__(self, Z, K, phi):
        self.z = tuple(Z)
        self.bigK = tuple(K)
        self.k = (K[0]*exp(Z[0]), K[1]*exp(Z[1]))
        self.phi = phi
        # The followint sets Ecoefficient and myC1C2
        calculate_coefficients(self)

    def gHat(self, a, b, c1, d1, c2, d2, s):
        z1, z2 = self.z
        phi = self.phi

        # Note: CheckList.m also uses [sigma, tau, q] but that function is not used anywhere.
        sigma = a/s**3 + b/s**2 -(a/2 + b + c1*exp(-z1) + c2*exp(-z2) )/s + (c1+d1)/(z1+s)+(c2+d2)/(z2+s)
        tau = a*(1/s**3 + 1/s**2) + b/s**2 -( z1*c1*exp(-z1)/(z1+s) + z2*c2*exp(-z2)/(z2+s) ) / s
        q = sigma - tau*exp(-s)

        ret = ( (a*(s+1) + b*s)/s**2 - z1*c1*exp(-z1)/(z1+s)-z2*c2*exp(-z2)/(z2+s) )*exp(-s) / (1-12*phi*q)
        return ret

    def CxCoef(self, d1, d2):
        z1, z2 = self.z
        k1, k2 = self.k

        phi = self.phi

        a, b, c1, c2 = self.ABC12(d1, d2)
        a00 = a**2
        b00 = -12*phi*((a+b)**2/2+a*( c1*exp(-z1)+c2*exp(-z2) ) )

        # Note: CxCoef.m had k(1)*exp(-Z(1))*exp(z(1))
        v1 = 24*phi*k1*self.gHat(a, b, c1, d1, c2, d2, z1)
        v2 = 24*phi*k2*self.gHat(a, b, c1, d1, c2, d2, z2)

        # TODO return a named tuple
        return [a00,b00,v1,v2,a,b,c1,d1,c2,d2]

    def Cr(self, a00, b00, v1, v2, r):
        phi = self.phi
        bigK1, bigK2 = self.bigK
        z1, z2 = self.z
        k1, k2 = self.k

        index = r <= 1
        rp = r[index]

        eCr = np.empty_like(r, dtype='d')
        eCr[index] = -(a00*rp+b00*rp**2 + 1/2 * phi * a00* rp**4 + v1/z1 *(1-exp(-z1*rp))+v2/z2*(1-exp(-z2*rp)) +
            v1**2/(2*k1*z1**2) *(cosh(z1*rp)-1)+ v2**2/(2*k2*z2**2) * (cosh(z2*rp)-1))

        index = r > 1
        rp = r[index]
        eCr[index] = (k1*exp(-z1*rp) + k2*exp(-z2*rp))

        eCr /= r
        return eCr

    def Ck(self, a00, b00, v1, v2, kk):
        bigK1, bigK2 = self.bigK
        z1, z2 = self.z
        phi = self.phi

        rou = phi*6/pi
        original = kk
        index = kk > 0.0001 - 1e-8
        kk = kk[index]

        def multiTerm(v, z, bigK, kk):
            return (v/z *( 1/kk**2 * (1-cos(kk)) - 1/(z**2 + kk**2)
    *(1-exp(-z)*(z/kk * sin(kk) + cos(kk)))) + v**2 / (4*bigK*z**2*exp(z)) *( 1/(z**2 + kk**2)
    *(2-exp(z)*(-z/kk * sin(kk) + cos(kk)) - exp(-z)*(z/kk * sin(kk) + cos(kk)))
    -2/kk**2 *(1-cos(kk))) - bigK/(z**2+kk**2)*(z/kk * sin(kk) + cos(kk)))

        eCk = np.empty_like(original, dtype='d')
        eCk[index] = (-24*phi*( a00*(-kk*cos(kk)+sin(kk))/kk**3 + b00/kk**4 *( -kk**2 * cos(kk) + 2*kk*sin(kk)
    +2*(cos(kk)-1)) +1/2 *phi*a00/kk**6 *( -kk**4*cos(kk) + 4*kk**3*sin(kk) + 12*kk**2*cos(kk)
    -24*kk*sin(kk)-24*(cos(kk)-1)) + multiTerm(v1,z1,bigK1,kk) + multiTerm(v2,z2,bigK2,kk)))

        index = original < 0.0001
        eCk[index] = (1-a00)
        eCk /= rou
        return eCk

    def hk(self, Ck):
        phi = self.phi
        rou = phi*6/pi
        hk = Ck / (1 - rou*Ck)
        return hk

    def Sk(self, hk):
        phi = self.phi
        rou = phi*6/pi
        return 1 + rou*hk

    def v(self, i):
        phi = self.phi
        z = self.z[i]
        ret = -12*phi*(1+2*phi)*(1-exp(-z)*(1+z))/z/(1-2*phi+phi**2) \
            +8*(1/z**2 - exp(-z)*(1/2 + (1+z)/z**2)) \
            +(8*(1+2*phi)*(-1+4*phi)*(1/z**2 - exp(-z)*(1/2+(1+z)/z**2)))/(1-2*phi+phi**2)
        return ret

    def w(self, i):
        phi = self.phi
        z = self.z[i]
        ret = 8/z**2 + 8*(1 + 2*phi)*(-1 + 4*phi)/(1-2*phi+phi**2)/z**2 - 12*phi*(1+2*phi)/(1-2*phi+phi**2)/z
        return ret

    def x(self, i):
        phi = self.phi
        z = self.z[i]
        ret=18*phi**2*(1-exp(-z)*(1+z))/z - \
            (12*phi*(-1+4*phi)*(1/z**2-exp(-z)*(1/2+(1+z)/z**2)))
        ret=ret/(1-2*phi+phi**2)
        return ret

    def y(self, i):
        phi = self.phi
        z = self.z[i]
        ret= -12*phi*(-1+4*phi)/(1-2*phi+phi**2)/z**2 + \
            18*phi**2/(1-2*phi+phi**2)/z
        return ret


def calculate_coefficients(self):

    z1, z2 = self.z
    k1, k2 = self.k
    phi = self.phi

    v1, v2 = self.v(0), self.v(1)
    w1, w2 = self.w(0), self.w(1)
    x1, x2 = self.x(0), self.x(1)
    y1, y2 = self.y(0), self.y(1)


    Ccd1_11 = -6*phi -6*exp(-2*z1)*phi+ 12*exp(-z1)*phi + 6*phi*v1 + 12*phi*x1 - 12*phi*v1/z1**2 + 12*exp(-z1)*phi*v1/z1**2 \
        +12*exp(-z1)*phi*v1/z1 - 12*phi*x1/z1 + 12*exp(-z1)*phi*x1/z1
    Ccd1_21 = 12*exp(-z2)*phi+6*phi*v2+12*phi*x2-12*phi*v2/z1**2+12*exp(-z1)*phi*v2/z1**2+12*exp(-z1)*phi*v2/z1-12*phi*x2/z1 \
        +12*exp(-z1)*phi*x2/z1-12*phi*z1/(z1+z2)-12*exp(-z2)*exp(-z1)*phi*z2/(z1+z2)
    Ccd2_12 = 12*exp(-z1)*phi+6*phi*v1+12*phi*x1-12*phi*v1/z2**2+12*exp(-z2)*phi*v1/z2**2+12*exp(-z2)*phi*v1/z2 \
        -12*phi*x1/z2+12*exp(-z2)*phi*x1/z2-12*phi*z2/(z1+z2)-12*exp(-z1-z2)*phi*z1/(z1+z2)
    Ccd2_22 = -6*phi+12*exp(-z2)*phi+6*phi*v2+12*phi*x2-12*phi*v2/z2**2+12*exp(-z2)*phi*v2/z2**2+12*exp(-z2)*phi*v2/z2 \
        -12*phi*x2/z2+12*exp(-z2)*phi*x2/z2-6*exp(-z2)*exp(-z2)*phi
    # TODO: "Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22" is a frequently used subexpression

    Cdd1_11 = -6*phi + 6*phi*w1 + 12*phi*y1-12*phi*w1/z1**2 + \
        12*exp(-z1)*phi*w1/z1**2 + 12*exp(-z1)*phi*w1/z1-12*phi*y1/z1+12*exp(-z1)*phi*y1/z1
    Cdd1_12 = 6*phi*w2+12*phi*y2-12*phi*w2/z1**2+12*exp(-z1)*phi*w2/z1**2+12*exp(-z1)*phi*w2/z1-12*phi*y2/z1+12*exp(-z1)*phi*y2/z1 \
        -12*phi*z1/(z1+z2)
    Cdd2_12 = 6*phi*w1+12*phi*y1-12*phi*w1/z2**2+12*exp(-z2)*phi*w1/z2**2+12*exp(-z2)*phi*w1/z2-12*phi*y1/z2 \
        +12*exp(-z2)*phi*y1/z2-12*phi*z2/(z1+z2)
    Cdd2_22 = -6*phi+6*phi*w2+12*phi*y2-12*phi*w2/z2**2+12*exp(-z2)*phi*w2/z2**2+12*exp(-z2)*phi*w2/z2-12*phi*y2/z2 \
        +12*exp(-z2)*phi*y2/z2

    a0 = (1 + 2*phi)/(1 - 2*phi+phi**2)
    b0 = -3*phi/2/(1-2*phi+phi**2)

    Cd1_1 = 6*a0*phi+12*b0*phi-12*a0*phi/z1**2+12*a0*exp(-z1)*phi/z1**2-12*b0*phi/z1+12*a0*exp(-z1)*phi/z1+12*b0*exp(-z1)*phi/z1+z1
    Cd2_2 = 6*a0*phi+12*b0*phi-12*a0*phi/z2**2+12*a0*exp(-z2)*phi/z2**2-12*b0*phi/z2+12*a0*exp(-z2)*phi/z2+12*b0*exp(-z2)*phi/z2+z2


    # Capture the intermediate variables for future calls to myC1C2
    def ABC12(d1, d2):
        A = ((aFNumd01*d2 + aFNumd21*d1**2*d2 + aFNumd10*d1 + aFNumd11*d1*d2 +
            aFNumd12*d1*d2**2)/(abFDend11*d1*d2))
        B = ((bNumd01*d2 + bFNumd21*d1**2*d2 + bFNumd10*d1 + bFNumd11*d1*d2 +
             bFNumd12*d1*d2**2)/(abFDend11 *d1*d2))
        C1 = ((-((Ccd2_22*Cd1_1*d1*d2 - Ccd1_21*Cd2_2*d1*d2 +
              Ccd2_22*Cdd1_11*d1**2*d2 -
              Ccd1_21*Cdd2_12*d1**2*d2 +
              Ccd2_22*Cdd1_12*d1*d2**2 -
              Ccd1_21*Cdd2_22*d1*d2**2 -
              Ccd2_22*d2*k1 + Ccd1_21*d1*k2)))/
        ((d1*((
              (-Ccd1_21)*Ccd2_12*d2 +
                Ccd1_11*Ccd2_22*d2)))))
        C2 = ((-((Ccd2_12*d2*
                (((-Cd1_1)*d1 - Cdd1_11*d1**2 -
                    Cdd1_12*d1*d2 + k1)) -
              Ccd1_11*d1*
                (((-Cd2_2)*d2 - Cdd2_12*d1*d2 -
                    Cdd2_22*d2**2 + k2))))) /
        (((-Ccd1_21)*Ccd2_12*d1*d2 +
            Ccd1_11*Ccd2_22*d1*d2)))
        return A, B, C1, C2
    self.ABC12 = ABC12

    aFNumd01 = -Ccd2_22*k1*v1+Ccd2_12*k1*v2
    aFNumd10 = Ccd1_21*k2*v1-Ccd1_11*k2*v2
    aFNumd11 = a0*Ccd1_21*Ccd2_12-a0*Ccd1_11*Ccd2_22+Ccd2_22*Cd1_1*v1-Ccd1_21*Cd2_2*v1-Ccd2_12*Cd1_1*v2+Ccd1_11*Cd2_2*v2
    aFNumd12 = Ccd2_22* Cdd1_12 * v1 - Ccd1_21* Cdd2_22* v1 - \
        Ccd2_12* Cdd1_12* v2 + Ccd1_11* Cdd2_22* v2 + \
        Ccd1_21* Ccd2_12* w2 - Ccd1_11* Ccd2_22* w2
    aFNumd21 = Ccd2_22*Cdd1_11*v1 - Ccd1_21*Cdd2_12*v1 - Ccd2_12*Cdd1_11*v2 + Ccd1_11*Cdd2_12*v2+Ccd1_21*Ccd2_12*w1 - Ccd1_11 * Ccd2_22*w1

    bFNumd10 = Ccd1_21*k2*x1 - Ccd1_11*k2*x2
    bFNumd11 = b0*Ccd1_21*Ccd2_12 - \
        b0*Ccd1_11*Ccd2_22 + Ccd2_22*Cd1_1*x1 - Ccd1_21*Cd2_2*x1 - Ccd2_12*Cd1_1*x2 + Ccd1_11*Cd2_2*x2
    bFNumd12 = Ccd2_22*Cdd1_12*x1 - Ccd1_21*Cdd2_22*x1 - \
        Ccd2_12*Cdd1_12*x2 + Ccd1_11 * Cdd2_22*x2 + Ccd1_21*Ccd2_12*y2 - Ccd1_11*Ccd2_22*y2
    bFNumd21 = Ccd2_22*Cdd1_11*x1 - Ccd1_21*Cdd2_12*x1 - \
        Ccd2_12*Cdd1_11*x2 + Ccd1_11*Cdd2_12*x2 + Ccd1_21*Ccd2_12*y1 - Ccd1_11*Ccd2_22*y1

    bNumd01 = -Ccd2_22*k1*x1+Ccd2_12*k1*x2
    abFDend11 = Ccd1_21*Ccd2_12-Ccd1_11*Ccd2_22

    c1F01 = Ccd2_22*k1
    c1F10 = -Ccd1_21*k2
    c1F11 = -Ccd2_22*Cd1_1 + Ccd1_21*Cd2_2
    c1F12 = -Ccd2_22*Cdd1_12 + Ccd1_21*Cdd2_22
    c1F21 = -Ccd2_22*Cdd1_11+Ccd1_21*Cdd2_12
    c2F01 = -Ccd2_12*k1
    c2F10 = Ccd1_11*k2
    c2F11 = Ccd2_12*Cd1_1 - Ccd1_11*Cd2_2
    c2F12 = Ccd2_12*Cdd1_12 - Ccd1_11 * Cdd2_22
    c2F21 = Ccd2_12*Cdd1_11 - Ccd1_11 * Cdd2_12
    c12FD = -Ccd1_21*Ccd2_12 + Ccd1_11*Ccd2_22

    def sigma_d01(s):
        return (aFNumd01/(abFDend11*s**3) + bNumd01/(abFDend11*s**2) -
                aFNumd01/(2*abFDend11*s) -
                bNumd01/(abFDend11*s) -
                (exp(-z2)*Ccd2_12*k1)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) +
                (exp(-z1)*Ccd2_22*k1)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) -
                (Ccd2_22*k1)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) +
                (Ccd2_12*k1)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z2)))

    def sigma_d10(s):
        return (aFNumd10/(abFDend11*s**3) + bFNumd10/(abFDend11*s**2) -
                aFNumd10/(2*abFDend11*s) -
                bFNumd10/(abFDend11*s) +
                (exp(-z2)*Ccd1_11*k2)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) -
                (exp(-z1)*Ccd1_21*k2)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) +
                (Ccd1_21*k2)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) -
                (Ccd1_11*k2)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z2)))

    def sigma_d11(s):
        return (aFNumd11/(abFDend11*s**3) + bFNumd11/(abFDend11*s**2) -
                aFNumd11/(2*abFDend11*s) -
                bFNumd11/(abFDend11*s) +
                (exp(-z2)*Ccd2_12*Cd1_1)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) -
                (exp(-z1)*Ccd2_22*Cd1_1)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) -
                (exp(-z2)*Ccd1_11*Cd2_2)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) +
                (exp(-z1)*Ccd1_21*Cd2_2)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) +
                (Ccd2_22*Cd1_1)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) -
                (Ccd1_21*Cd2_2)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) -
                (Ccd2_12*Cd1_1)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z2)) +
                (Ccd1_11*Cd2_2)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z2)))

    def sigma_d12(s):
        return (aFNumd12/(abFDend11*s**3) +
                bFNumd12/(abFDend11*s**2) -
                aFNumd12/(2*abFDend11*s) -
                bFNumd12/(abFDend11*s) +
                (exp(-z2)*Ccd2_12*Cdd1_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) -
                (exp(-z1)*Ccd2_22*Cdd1_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) -
                (exp(-z2)*Ccd1_11*Cdd2_22)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) +
                (exp(-z1)*Ccd1_21*Cdd2_22)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) +
                (Ccd2_22*Cdd1_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) -
                (Ccd1_21*Cdd2_22)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) +
                1/(s + z2) -
                (Ccd2_12*Cdd1_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z2)) +
                (Ccd1_11*Cdd2_22)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z2)))

    def sigma_d21(s):
        return (aFNumd21/(abFDend11*s**3) + bFNumd21/(abFDend11*s**2) -
                aFNumd21/(2*abFDend11*s) -
                bFNumd21/(abFDend11*s) +
                (exp(-z2)*Ccd2_12*Cdd1_11)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) -
                (exp(-z1)*Ccd2_22*Cdd1_11)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) -
                (exp(-z2)*Ccd1_11*Cdd2_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) +
                (exp(-z1)*Ccd1_21*Cdd2_12)/(s*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)) +
                (Ccd1_21*Ccd2_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) -
                (Ccd1_11*Ccd2_22)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) +
                (Ccd2_22*Cdd1_11)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) -
                (Ccd1_21*Cdd2_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z1)) -
                (Ccd2_12*Cdd1_11)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z2)) +
                (Ccd1_11*Cdd2_12)/((Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*(s + z2)))

    def tau_d01(s):
        return (aFNumd01/(abFDend11*s**3) + aFNumd01/(abFDend11*s**2) +
                bNumd01/(abFDend11*s**2) +
                (exp(-z1)*Ccd2_22*k1*z1)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z1)) -
                (exp(-z2)*Ccd2_12*k1*z2)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z2)))

    def tau_d10(s):
        return (aFNumd10/(abFDend11*s**3) + aFNumd10/(abFDend11*s**2) +
                bFNumd10/(abFDend11*s**2) -
                (exp(-z1)*Ccd1_21*k2*z1)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z1)) +
                (exp(-z2)*Ccd1_11*k2*z2)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z2)))

    def tau_d11(s):
        return (aFNumd11/(abFDend11*s**3) + aFNumd11/(abFDend11*s**2) +
                bFNumd11/(abFDend11*s**2) -
                (exp(-z1)*Ccd2_22*Cd1_1*z1)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z1)) +
                (exp(-z1)*Ccd1_21*Cd2_2*z1)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z1)) +
                (exp(-z2)*Ccd2_12*Cd1_1*z2)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z2)) -
                (exp(-z2)*Ccd1_11*Cd2_2*z2)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z2)))

    def tau_d12(s):
        return (aFNumd12/(abFDend11*s**3) + aFNumd12/(abFDend11*s**2) +
                bFNumd12/(abFDend11*s**2) -
                (exp(-z1)*Ccd2_22*Cdd1_12*z1)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z1)) +
                (exp(-z1)*Ccd1_21*Cdd2_22*z1)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z1)) +
                (exp(-z2)*Ccd2_12*Cdd1_12*z2)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z2)) -
                (exp(-z2)*Ccd1_11*Cdd2_22*z2)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z2)))

    def tau_d21(s):
        return (aFNumd21/(abFDend11*s**3) + aFNumd21/(abFDend11*s**2) +
                bFNumd21/(abFDend11*s**2) -
                (exp(-z1)*Ccd2_22*Cdd1_11*z1)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z1)) +
                (exp(-z1)*Ccd1_21*Cdd2_12*z1)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z1)) +
                (exp(-z2)*Ccd2_12*Cdd1_11*z2)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z2)) -
                (exp(-z2)*Ccd1_11*Cdd2_12*z2)/(s*
                    (Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)*
                    (s + z2)))

    E1d02 = 12*c1F01*phi*sigma_d01(z1) - 12*c1F01*exp(-z1)*phi*tau_d01(z1)

    E1d11 = (((12*c1F10*phi*sigma_d01(z1) + \
          12*c1F01*phi*sigma_d10(z1) - \
          12*c1F10*exp(-z1)*phi*tau_d01(z1) - \
          12*c1F01*exp(-z1)*phi*tau_d10(z1))))
    E1d12 = (-c1F01 + 12*c1F11*phi*sigma_d01(z1) +
             12*c1F01*phi*sigma_d11(z1) -
             12*c1F11*exp(-z1)*phi*tau_d01(z1) -
             12*c1F01*exp(-z1)*phi*tau_d11(z1))
    E1d13 = (12*c1F12*phi*sigma_d01(z1) +
             12*c1F01*phi*sigma_d12(z1) -
             12*c1F12*exp(-z1)*phi*tau_d01(z1) -
             12*c1F01*exp(-z1)*phi*tau_d12(z1))

    E1d20 = 12*c1F10*phi*sigma_d10(z1) - 12*c1F10*exp(-z1)*phi*tau_d10(z1)
    E1d21 = (-c1F10 + 12*c1F11*phi*sigma_d10(z1) +
             12*c1F10*phi*sigma_d11(z1) -
             12*c1F11*exp(-z1)*phi*tau_d10(z1) -
             12*c1F10*exp(-z1)*phi*tau_d11(z1))
    E1d22 = (-c1F11 + 12*c12FD*phi*sigma_d01(z1) +
             12*c1F21*phi*sigma_d01(z1) +
             12*c1F12*phi*sigma_d10(z1) +
             12*c1F11*phi*sigma_d11(z1) +
             12*c1F10*phi*sigma_d12(z1) +
             12*c1F01*phi*sigma_d21(z1) -
             12*c1F21*exp(-z1)*phi*tau_d01(z1) -
             12*c1F12*exp(-z1)*phi*tau_d10(z1) -
             12*c1F11*exp(-z1)*phi*tau_d11(z1) -
             12*c1F10*exp(-z1)*phi*tau_d12(z1) -
             12*c1F01*exp(-z1)*phi*tau_d21(z1))
    E1d23 = (-c1F12 + 12*c1F12*phi*sigma_d11(z1) +
             12*c1F11*phi*sigma_d12(z1) -
             12*c1F12*exp(-z1)*phi*tau_d11(z1) -
             12*c1F11*exp(-z1)*phi*tau_d12(z1))
    E1d24 = 12*c1F12*phi*sigma_d12(z1) - 12*c1F12*exp(-z1)*phi*tau_d12(z1)

    E1d31 = (12*c12FD*phi*sigma_d10(z1) +
             12*c1F21*phi*sigma_d10(z1) +
             12*c1F10*phi*sigma_d21(z1) -
             12*c1F21*exp(-z1)*phi*tau_d10(z1) -
             12*c1F10*exp(-z1)*phi*tau_d21(z1))
    E1d32 = (-c12FD - c1F21 + 12*c12FD*phi*sigma_d11(z1) +
             12*c1F21*phi*sigma_d11(z1) +
             12*c1F11*phi*sigma_d21(z1) -
             12*c1F21*exp(-z1)*phi*tau_d11(z1) -
             12*c1F11*exp(-z1)*phi*tau_d21(z1))
    E1d33 = (12*c12FD*phi*sigma_d12(z1) +
             12*c1F21*phi*sigma_d12(z1) +
             12*c1F12*phi*sigma_d21(z1) -
             12*c1F21*exp(-z1)*phi*tau_d12(z1) -
             12*c1F12*exp(-z1)*phi*tau_d21(z1))

    E1d42 = (12*c12FD*phi*sigma_d21(z1) +
             12*c1F21*phi*sigma_d21(z1) -
             12*c1F21*exp(-z1)*phi*tau_d21(z1))

    E2d02 = 12*c2F01*exp(z2)*phi*sigma_d01(z2)-12*c2F01*phi*tau_d01(z2)

    E2d11 = (12*c2F10*exp(z2)*phi*sigma_d01(z2) +
             12*c2F01*exp(z2)*phi*sigma_d10(z2) -
             12*c2F10*phi*tau_d01(z2) -
             12*c2F01*phi*tau_d10(z2))
    E2d12 = (-c2F01*exp(z2) +
             12*c2F11*exp(z2)*phi*sigma_d01(z2) +
             12*c2F01*exp(z2)*phi*sigma_d11(z2) -
             12*c2F11*phi*tau_d01(z2) -
             12*c2F01*phi*tau_d11(z2))
    E2d13 = (12*c12FD*exp(z2)*phi*sigma_d01(z2) +
             12*c2F12*exp(z2)*phi*sigma_d01(z2) +
             12*c2F01*exp(z2)*phi*sigma_d12(z2) -
             12*c2F12*phi*tau_d01(z2) -
             12*c2F01*phi*tau_d12(z2))

    E2d20 = 12*c2F10*exp(z2)*phi*sigma_d10(z2) - 12*c2F10*phi*tau_d10(z2)
    E2d21 = (-c2F10*exp(z2) +
             12*c2F11*exp(z2)*phi*sigma_d10(z2) +
             12*c2F10*exp(z2)*phi*sigma_d11(z2) -
             12*c2F11*phi*tau_d10(z2) -
             12*c2F10*phi*tau_d11(z2))
    E2d22 = (-c2F11*exp(z2) +
             12*c2F21*exp(z2)*phi*sigma_d01(z2) +
             12*c12FD*exp(z2)*phi*sigma_d10(z2) +
             12*c2F12*exp(z2)*phi*sigma_d10(z2) +
             12*c2F11*exp(z2)*phi*sigma_d11(z2) +
             12*c2F10*exp(z2)*phi*sigma_d12(z2) +
             12*c2F01*exp(z2)*phi*sigma_d21(z2) -
             12*c2F21*phi*tau_d01(z2) -
             12*c2F12*phi*tau_d10(z2) -
             12*c2F11*phi*tau_d11(z2) -
             12*c2F10*phi*tau_d12(z2) -
             12*c2F01*phi*tau_d21(z2))
    E2d23 = (-c12FD*exp(z2) - c2F12*exp(z2) +
             12*c12FD*exp(z2)*phi*sigma_d11(z2) +
             12*c2F12*exp(z2)*phi*sigma_d11(z2) +
             12*c2F11*exp(z2)*phi*sigma_d12(z2) -
             12*c2F12*phi*tau_d11(z2) -
             12*c2F11*phi*tau_d12(z2))
    E2d24 = (12*c12FD*exp(z2)*phi*sigma_d12(z2) +
             12*c2F12*exp(z2)*phi*sigma_d12(z2) -
             12*c2F12*phi*tau_d12(z2))

    E2d31 = (12*c2F21*exp(z2)*phi*sigma_d10(z2) +
             12*c2F10*exp(z2)*phi*sigma_d21(z2) -
             12*c2F21*phi*tau_d10(z2) -
             12*c2F10*phi*tau_d21(z2))
    E2d32 = (-c2F21*exp(z2) +
             12*c2F21*exp(z2)*phi*sigma_d11(z2) +
             12*c2F11*exp(z2)*phi*sigma_d21(z2) -
             12*c2F21*phi*tau_d11(z2) -
             12*c2F11*phi*tau_d21(z2))
    E2d33 = (12*c2F21*exp(z2)*phi*sigma_d12(z2) +
             12*c12FD*exp(z2)*phi*sigma_d21(z2) +
             12*c2F12*exp(z2)*phi*sigma_d21(z2) -
             12*c2F21*phi*tau_d12(z2) -
             12*c2F12*phi*tau_d21(z2))

    E2d42 = (12*c2F21*exp(z2)*phi*sigma_d21(z2) -
             12*c2F21*phi*tau_d21(z2))


    d1Factor, d2Factor = self.d1Factor, self.d2Factor

    # Assign values to Ecoefficient array
    Ecoefficient = np.empty(26, 'd')
    Ecoefficient[0] = E1d20 * 0
    Ecoefficient[1] = E1d11
    Ecoefficient[2] = E1d21 * 0
    Ecoefficient[3] = E1d31
    Ecoefficient[4] = E1d02
    Ecoefficient[5] = E1d12
    Ecoefficient[22] = E1d22
    Ecoefficient[23] = E1d32
    Ecoefficient[6] = E1d42
    Ecoefficient[7] = E1d13
    Ecoefficient[8] = E1d23 * 0
    Ecoefficient[9] = E1d33
    Ecoefficient[10] = E1d24 * 0

    Ecoefficient[11] = E2d20
    Ecoefficient[12] = E2d11
    Ecoefficient[13] = E2d21
    Ecoefficient[14] = E2d31
    Ecoefficient[15] = E2d02 * 0
    Ecoefficient[16] = E2d12 * 0
    Ecoefficient[24] = E2d22
    Ecoefficient[25] = E2d32 * 0
    Ecoefficient[17] = E2d42 * 0
    Ecoefficient[18] = E2d13
    Ecoefficient[19] = E2d23
    Ecoefficient[20] = E2d33
    Ecoefficient[21] = E2d24

    # Apply scaling factors
    Ecoefficient[0] = Ecoefficient[0] * d1Factor**2 * 1
    Ecoefficient[1] = Ecoefficient[1] * d1Factor * d2Factor
    Ecoefficient[2] = Ecoefficient[2] * d1Factor**2 * d2Factor
    Ecoefficient[3] = Ecoefficient[3] * d1Factor**3 * d2Factor
    Ecoefficient[4] = Ecoefficient[4] * 1 * d2Factor**2
    Ecoefficient[5] = Ecoefficient[5] * d1Factor**1 * d2Factor**2
    Ecoefficient[22] = Ecoefficient[22] * d1Factor**2 * d2Factor**2
    Ecoefficient[23] = Ecoefficient[23] * d1Factor**3 * d2Factor**2
    Ecoefficient[6] = Ecoefficient[6] * d1Factor**4 * d2Factor**2
    Ecoefficient[7] = Ecoefficient[7] * d1Factor**1 * d2Factor**3
    Ecoefficient[8] = Ecoefficient[8] * d1Factor**2 * d2Factor**3
    Ecoefficient[9] = Ecoefficient[9] * d1Factor**3 * d2Factor**3
    Ecoefficient[10] = Ecoefficient[10] * d1Factor**2 * d2Factor**4

    Ecoefficient[11] = Ecoefficient[11] * d1Factor**2 * 1
    Ecoefficient[12] = Ecoefficient[12] * d1Factor**1 * d2Factor**1
    Ecoefficient[13] = Ecoefficient[13] * d1Factor**2 * d2Factor**1
    Ecoefficient[14] = Ecoefficient[14] * d1Factor**3 * d2Factor**1
    Ecoefficient[15] = Ecoefficient[15] * 1 * d2Factor**2
    Ecoefficient[16] = Ecoefficient[16] * d1Factor**1 * d2Factor**2
    Ecoefficient[24] = Ecoefficient[24] * d1Factor**2 * d2Factor**2
    Ecoefficient[25] = Ecoefficient[25] * d1Factor**3 * d2Factor**2
    Ecoefficient[17] = Ecoefficient[17] * d1Factor**4 * d2Factor**2
    Ecoefficient[18] = Ecoefficient[18] * d1Factor**1 * d2Factor**3
    Ecoefficient[19] = Ecoefficient[19] * d1Factor**2 * d2Factor**3
    Ecoefficient[20] = Ecoefficient[20] * d1Factor**3 * d2Factor**3
    Ecoefficient[21] = Ecoefficient[21] * d1Factor**2 * d2Factor**4

    self.Ecoefficient = Ecoefficient

    # The following are used in d_coeI_J functions, but those aren't used. Suppress for now.
    """
    def aDend1(d2): return d2*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)
    def aNumd0(d2): return d2*(-Ccd2_22*k1*v1+Ccd2_12*k1*v2)
    def aNumd1(d2): return (a0*Ccd1_21*Ccd2_12*d2 - a0*Ccd1_11*Ccd2_22*d2 -
        Ccd1_21*Cdd2_12*v1 + Ccd2_22*Cd1_1*d2*v1 -
        Ccd1_21*Cd2_2*d2*v1 +
        Ccd2_22*Cdd1_12*d2**2*v1 -
        Ccd1_21*Cdd2_22*d2**2*v1 + Ccd1_21*k2*v1 +
        Ccd1_11*Cdd2_12*v2 - Ccd2_12*Cd1_1*d2*v2 +
        Ccd1_11*Cd2_2*d2*v2 -
        Ccd2_12*Cdd1_12*d2**2*v2 +
        Ccd1_11*Cdd2_22*d2**2*v2 - Ccd1_11*k2*v2 +
        Ccd1_21*Ccd2_12*d2**2*w2 -
        Ccd1_11*Ccd2_22*d2**2*w2)
    def aNumd2(d2): return d2*(Ccd2_22*Cdd1_11*v1 - Ccd2_12*Cdd1_11*v2 + Ccd1_21*Ccd2_12*w1 - Ccd1_11*Ccd2_22*w1)

    def bDend1(d2): return d2*(Ccd1_21*Ccd2_12-Ccd1_11*Ccd2_22)
    def bNumd0(d2): return d2*(-Ccd2_22*k1*x1 + Ccd2_12*k1*x2)
    def bNumd1(d2): return (b0*Ccd1_21*Ccd2_12*d2 - b0*Ccd1_11*Ccd2_22*d2 -
        Ccd1_21*Cdd2_12*x1 + Ccd2_22*Cd1_1*d2*x1 -
        Ccd1_21*Cd2_2*d2*x1 +
        Ccd2_22*Cdd1_12*d2**2*x1 -
        Ccd1_21*Cdd2_22*d2**2*x1 + Ccd1_21*k2*x1 +
        Ccd1_11*Cdd2_12*x2 - Ccd2_12*Cd1_1*d2*x2 +
        Ccd1_11*Cd2_2*d2*x2 -
        Ccd2_12*Cdd1_12*d2**2*x2 +
        Ccd1_11*Cdd2_22*d2**2*x2 - Ccd1_11*k2*x2 +
        Ccd1_21*Ccd2_12*d2**2*y2 -
        Ccd1_11*Ccd2_22*d2**2*y2)
    def bNumd2(d2): return d2*(Ccd2_22*Cdd1_11*x1 - Ccd2_12*Cdd1_11*x2 + Ccd1_21*Ccd2_12*y1 - Ccd1_11*Ccd2_22*y1)

    def cDend1_1(d2): return d2*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)
    def cDend1_2(d2): return d2*(Ccd1_21*Ccd2_12 - Ccd1_11*Ccd2_22)
    def cNumd0_1(d2): return d2*(-Ccd2_22*d2*k1)
    def cNumd0_2(d2): return d2*(Ccd2_12*k1)
    def cNumd1_1(d2): return d2*(Ccd2_22*Cd1_1 - Ccd1_21*Cd2_2 + Ccd2_22*Cdd1_12*d2 - Ccd1_21*Cdd2_22*d2 + Ccd1_21*k2)
    def cNumd1_2(d2): return d2*(-Ccd2_12*Cd1_1 + Ccd1_11*Cd2_2 - Ccd2_12*Cdd1_12*d2 + Ccd1_11*Cdd2_22*d2 - Ccd1_11*k2)
    def cNumd2_1(d2): return d2*(Ccd2_22*Cdd1_11 - Ccd1_21*Cdd2_12)
    def cNumd2_2(d2): return d2*(-Ccd2_12*Cdd1_11 + Ccd1_11*Cdd2_12)
    """

