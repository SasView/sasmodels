from __future__ import division, print_function

import numpy as np
from numpy import sqrt, pi, tan, cos, sin, log, exp, arctan2 as atan2, sign, radians, degrees

# mpmath version of special functions
_ = """
import mpmath as mp
sn, cn, dn = [mp.ellipfun(v) for v in 'sn', 'cn', 'dn']

def ellipticJi(u, v, m):
    z = u+v*1j
    return sn(z, m), cn(z, m), dn(z, m)

def ellipticJ(u, m):
    s, c, d = sn(u, m), cn(u, m), dn(u, m)
    phi = mp.asin(s)
    return s, c, d, phi

def ellipticFi(phi, psi, m):
    z = phi + psi*1j
    return mp.ellipf(z, m)

def ellipticF(phi, m):
    return mp.ellipf(phi, m)
"""

# scipy version of special functions
from scipy.special import ellipj as ellipticJ, ellipkinc as ellipticF
from numpy import sinh, sign, arctan as atan

def ellipticJi(u, v, m):
    scalar = np.isscalar(u) and np.isscalar(v) and np.isscalar(m)
    u, v, m = np.broadcast_arrays(u, v, m)
    result = np.empty_like([u,u,u], 'D')
    real = v==0
    imag = u==0
    mixed = ~(real|imag)
    result[:, real] = _ellipticJi_real(u[real], m[real])
    result[:, imag] = _ellipticJi_imag(v[imag], m[imag])
    result[:, mixed] = _ellipticJi(u[mixed], v[mixed], m[mixed])
    return result[0,:] if scalar else result

def _ellipticJi_real(u, m):
    sn, cn, dn, phi = ellipticJ(u, m)
    return sn, cn, dn

def _ellipticJi_imag(v, m):
    sn, cn, dn, phi = ellipticJ(v, 1-m)
    return 1j*sn/cn, 1/cn, dn/cn

def _ellipticJi(u, v, m):
    # Ignoring special cases for now.
    #   u=0: (1j*b[0]/b[1], 1/b[1], b[2]/b[1])
    #   v=0: (a[0], a[1], a[2])
    a = ellipticJ(u, m)
    b = ellipticJ(v, 1 - m)
    c = b[1]**2 + m * (a[0] * b[0])**2
    return [
        (a[0] * b[2] / c) + 1j*(a[1] * a[2] * b[0] * b[1] / c),
        (a[1] * b[1] / c) + 1j*(-a[0] * a[2] * b[0] * b[2] / c),
        (a[2] * b[1] * b[2] / c) + 1j*(-m * a[0] * a[1] * b[0] / c),
    ]

# calculate F(phi+ipsi|m).
# see Abramowitz and Stegun, 17.4.11.
def ellipticFi(phi, psi, m):
    if np.any(phi == 0):
        scalar = np.isscalar(phi) and np.isscalar(psi) and np.isscalar(m)
        phi, psi, m = np.broadcast_arrays(phi, psi, m)
        result = np.empty_like(phi, 'D')
        index = (phi == 0)
        result[index] = ellipticF(atan(sinh(abs(phi[index]))), 1-m[index]) * sign(psi[index])
        result[~index] = ellipticFi(phi[~index], psi[~index], m[~index])
        return result.reshape(1)[0] if scalar else result

    r = abs(phi)
    i = abs(psi)
    sinhpsi2 = sinh(i)**2
    cscphi2 = 1 / sin(r)**2
    cotphi2 = 1 / tan(r)**2
    b = -(cotphi2 + m * (sinhpsi2 * cscphi2) - 1 + m)
    c = (m - 1) * cotphi2
    cotlambda2 = (-b + sqrt(b * b - 4 * c)) / 2
    re = ellipticF(atan(1 / sqrt(cotlambda2)), m) * sign(phi)
    im = ellipticF(atan(sqrt((cotlambda2 / cotphi2 - 1) / m)), 1 - m) * sign(psi)
    return re + 1j*im

_ = """
// Returns [sn, cn, dn](u + iv|m).
export function ellipticJi(u, v, m) {
  var a, b, c;
  if (!u) {
    b = ellipticJ(v, 1 - m);
    return [
      [0, b[0] / b[1]],
      [1 / b[1], 0],
      [b[2] / b[1], 0]
    ];
  }
  a = ellipticJ(u, m);
  if (!v) return [[a[0], 0], [a[1], 0], [a[2], 0]];
  b = ellipticJ(v, 1 - m);
  c = b[1] * b[1] + m * a[0] * a[0] * b[0] * b[0];
  return [
    [a[0] * b[2] / c, a[1] * a[2] * b[0] * b[1] / c],
    [a[1] * b[1] / c, -a[0] * a[2] * b[0] * b[2] / c],
    [a[2] * b[1] * b[2] / c, -m * a[0] * a[1] * b[0] / c]
  ];
}

// Returns [sn, cn, dn, ph](u|m).
export function ellipticJ(u, m) {
  var ai, b, phi, t, twon;
  if (m < epsilon) {
    t = sin(u);
    b = cos(u);
    ai = m * (u - t * b) / 4;
    return [
      t - ai * b,
      b + ai * t,
      1 - m * t * t / 2,
      u - ai
    ];
  }
  if (m >= 1 - epsilon) {
    ai = (1 - m) / 4;
    b = cosh(u);
    t = tanh(u);
    phi = 1 / b;
    twon = b * sinh(u);
    return [
      t + ai * (twon - u) / (b * b),
      phi - ai * t * phi * (twon - u),
      phi + ai * t * phi * (twon + u),
      2 * atan(exp(u)) - halfPi + ai * (twon - u) / b
    ];
  }

  var a = [1, 0, 0, 0, 0, 0, 0, 0, 0],
      c = [sqrt(m), 0, 0, 0, 0, 0, 0, 0, 0],
      i = 0;
  b = sqrt(1 - m);
  twon = 1;

  while (abs(c[i] / a[i]) > epsilon && i < 8) {
    ai = a[i++];
    c[i] = (ai - b) / 2;
    a[i] = (ai + b) / 2;
    b = sqrt(ai * b);
    twon *= 2;
  }

  phi = twon * a[i] * u;
  do {
    t = c[i] * sin(b = phi) / a[i];
    phi = (asin(t) + phi) / 2;
  } while (--i);

  return [sin(phi), t = cos(phi), t / cos(phi - b), phi];
}

// calculate F(phi+ipsi|m).
// see Abramowitz and Stegun, 17.4.11.
export function ellipticFi(phi, psi, m) {
  var r = abs(phi),
      i = abs(psi),
      sinhpsi = sinh(i);
  if (r) {
    var cscphi = 1 / sin(r),
        cotphi2 = 1 / (tan(r) * tan(r)),
        b = -(cotphi2 + m * (sinhpsi * sinhpsi * cscphi * cscphi) - 1 + m),
        c = (m - 1) * cotphi2,
        cotlambda2 = (-b + sqrt(b * b - 4 * c)) / 2;
    return [
      ellipticF(atan(1 / sqrt(cotlambda2)), m) * sign(phi),
      ellipticF(atan(sqrt((cotlambda2 / cotphi2 - 1) / m)), 1 - m) * sign(psi)
    ];
  }
  return [
    0,
    ellipticF(atan(sinhpsi), 1 - m) * sign(psi)
  ];
}

// Calculate F(phi|m) where m = k^2 = sin(alpha)^2.
// See Abramowitz and Stegun, 17.6.7.
export function ellipticF(phi, m) {
  if (!m) return phi;
  if (m === 1) return log(tan(phi / 2 + quarterPi));
  var a = 1,
      b = sqrt(1 - m),
      c = sqrt(m);
  for (var i = 0; abs(c) > epsilon; i++) {
    if (phi % pi) {
      var dPhi = atan(b * tan(phi) / a);
      if (dPhi < 0) dPhi += pi;
      phi += dPhi + ~~(phi / pi) * pi;
    } else phi += phi;
    c = (a + b) / 2;
    b = sqrt(a * b);
    c = ((a = c) - b) / 2;
  }
  return phi / (pow(2, i) * a);
"""

def guyouComplexAtan(x, y):
    x2 = x * x
    y_1 = y + 1
    t = 1 - x2 - y * y
    return (
       0.5 * (sign(x+(x==0))*halfPi - atan2(t, 2 * x)),
       -0.25 * log(t * t + 4 * x2) +0.5 * log(y_1 * y_1 + x2)
    )

sqrt2 = sqrt(2)
halfPi = pi/2

# From d3-geo-projections
def guyou(lam, phi):
    x, y = np.asarray(lam), np.asarray(phi)
    xn, x = divmod(x+90, 180)
    yn, y = divmod(y+90, 180)
    xn, lam = xn*180, radians(x-90)
    yn, phi = yn*180, radians(y-90)

    k_ = (sqrt2 - 1) / (sqrt2 + 1)
    k = sqrt(1 - k_ * k_)
    K = ellipticF(halfPi, k * k)
    f = -1
    psi = log(tan(pi / 4 + abs(phi) / 2))
    r = exp(f * psi) / sqrt(k_)
    at = guyouComplexAtan(r * cos(f * lam), r * sin(f * lam))
    t = ellipticFi(at[0], at[1], k * k)
    #return [-t.imag, (1 if phi >=0 else -1)*(0.5 * K - t.real)]
    x, y = (-t.imag, sign(phi + (phi==0))*(0.5 * K - t.real))
    return degrees(x)+xn, degrees(y)+yn

#def guyouComplexDivide(a, b):
#    denominator = b[0] * b[0] + b[1] * b[1]
#    return [
#        (a[0] * b[0] + a[1] * b[1]) / denominator,
#        (a[1] * b[0] - a[0] * b[1]) / denominator
#    ]

def guyou_invert(x, y):
    x, y = np.asarray(x), np.asarray(y)
    xn, x = divmod(x+90, 180)
    yn, y = divmod(y+90, 180)
    xn, x = xn*180, radians(x-90)
    yn, y = yn*180, radians(y-90)

    k_ = (sqrt2 - 1) / (sqrt2 + 1)
    k = sqrt(1 - k_ * k_)
    K = ellipticF(halfPi, k * k)
    f = -1
    j = ellipticJi(0.5 * K - y, -x, k * k)
    #tn = guyouComplexDivide(j[0], j[1])
    #lam = atan2(tn[1], tn[0]) / f
    tn = j[0]/j[1]  # j[0], j[1] are complex
    lam = atan2(tn.imag, tn.real) / f
    phi = 2 * atan(exp(0.5 / f * log(k_ * tn.real**2 + k_ * tn.imag**2))) - halfPi
    return degrees(lam)+xn, degrees(phi)+yn

def plot_grid():
    import matplotlib.pyplot as plt
    from numpy import linspace
    lat_line = linspace(-90, 90, 100)[1:-1]
    long_line = linspace(-90, 90, 100)[1:-1]

    #scale = 1
    scale = 2
    plt.subplot(211)
    for lat in range(-80, 81, 10):
        x, y = guyou(scale*lat_line, scale*lat)
        plt.plot(x, y, 'g')

    for long in range(-80, 81, 10):
        x, y = guyou(scale*long, scale*long_line)
        plt.plot(x, y, 'b')
    plt.subplot(212)
    for lat in range(-80, 81, 10):
        x, y = guyou_invert(scale*lat_line, scale*lat)
        plt.plot(x, y, 'g')

    for long in range(-80, 81, 10):
        x, y = guyou_invert(scale*long, scale*long_line)
        plt.plot(x, y, 'b')

if __name__ == "__main__":
    plot_grid()
    import matplotlib.pyplot as plt; plt.show()
