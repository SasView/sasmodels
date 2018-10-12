# Copyright 2013-2016 Mike Bostock
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the author nor the names of contributors may be used to
#   endorse or promote products derived from this software without specific prior
#   written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# https://github.com/d3/d3-geo-projection
# commit fd2886555e46b35163c7b898d43c7d1bcebbba7c 2016-07-02
#
# 2017-11-01 Paul Kienzle
# * converted to python, with degrees rather than radians
"""
Convert between latitude-longitude and Guyou map coordinates.
"""

from __future__ import division, print_function

import numpy as np
from numpy import sqrt, pi, tan, cos, sin, sign, radians, degrees
from numpy import sinh, arctan as atan

# scipy version of special functions
from scipy.special import ellipj as ellipticJ, ellipkinc as ellipticF

_ = """
# mpmath version of special functions
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

def ellipticJi(u, v, m):
    scalar = np.isscalar(u) and np.isscalar(v) and np.isscalar(m)
    u, v, m = np.broadcast_arrays(u, v, m)
    result = np.empty_like([u, u, u], 'D')
    real = (v == 0)
    imag = (u == 0)
    mixed = ~(real|imag)
    result[:, real] = _ellipticJi_real(u[real], m[real])
    result[:, imag] = _ellipticJi_imag(v[imag], m[imag])
    result[:, mixed] = _ellipticJi(u[mixed], v[mixed], m[mixed])
    return result[0, :] if scalar else result

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
        result[index] = ellipticF(atan(sinh(abs(phi[index]))),
                                  1-m[index]) * sign(psi[index])
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
    im = ellipticF(atan(sqrt(np.maximum(0, (cotlambda2 / cotphi2 - 1) / m))),
                   1 - m) * sign(psi)
    return re + 1j*im

SQRT2 = sqrt(2)

# [PAK] renamed k_ => cos_u, k => sin_u, k*k => sinsq_u to avoid k,K confusion
# cos_u = 0.171572875253809902396622551580603842860656249246103853646...
# sinsq_u = 0.970562748477140585620264690516376942836062504523376878120...
# K = 3.165103454447431823666270142140819753058976299237578486994...
def guyou(lam, phi):
    """Transform from (latitude, longitude) to point (x, y)"""
    # [PAK] wrap into [-pi/2, pi/2] radians
    x, y = np.asarray(lam), np.asarray(phi)
    xn, x = divmod(x+90, 180)
    yn, y = divmod(y+90, 180)
    xn, lam = xn*180, radians(x-90)
    yn, phi = yn*180, radians(y-90)

    # Compute constant K
    cos_u = (SQRT2 - 1) / (SQRT2 + 1)
    sinsq_u = 1 - cos_u**2
    K = ellipticF(pi/2, sinsq_u)

    # [PAK] simplify expressions, using the fact that f = -1
    # Note: exp(f log(x)) => 1/x,  cos(f x) => cos(x), sin(f x) => -sin(x)
    r = 1/(tan(pi/4 + abs(phi)/2) * sqrt(cos_u))
    at = atan(r * (cos(lam) - 1j*sin(lam)))
    t = ellipticFi(at.real, at.imag, sinsq_u)
    x, y = (-t.imag, sign(phi + (phi == 0))*(0.5 * K - t.real))

    # [PAK] convert to degrees, and return to original tile
    return degrees(x)+xn, degrees(y)+yn

def guyou_invert(x, y):
    """Transform from point (x, y) on plot to (latitude, longitude)"""
    # [PAK] wrap into [-pi/2, pi/2] radians
    x, y = np.asarray(x), np.asarray(y)
    xn, x = divmod(x+90, 180)
    yn, y = divmod(y+90, 180)
    xn, x = xn*180, radians(x-90)
    yn, y = yn*180, radians(y-90)

    # compute constant K
    cos_u = (SQRT2 - 1) / (SQRT2 + 1)
    sinsq_u = 1 - cos_u**2
    K = ellipticF(pi/2, sinsq_u)

    # [PAK] simplify expressions, using the fact that f = -1
    j = ellipticJi(K/2 - y, -x, sinsq_u)
    tn = j[0]/j[1]  # j[0], j[1] are complex
    # Note: -atan2(im(x)/re(x)) => angle(x)
    lam = -np.angle(tn)
    # Note: exp(0.5/f log(a re(x)^2 + a im(x)^2)) => 1/(sqrt(a) |x|)
    phi = 2*atan(1/sqrt(cos_u)/abs(tn)) - pi/2

    # [PAK] convert to degrees, and return to original tile
    return degrees(lam)+xn, degrees(phi)+yn

def plot_grid():
    """Plot the latitude-longitude grid for Guyou transform"""
    import matplotlib.pyplot as plt
    from numpy import linspace
    lat_line = linspace(-90, 90, 400)
    long_line = linspace(-90, 90, 400)

    #scale = 1
    limit, step, scale = 90, 10, 2
    plt.subplot(211)
    for lat in range(-limit, limit+1, step):
        x, y = guyou(scale*lat_line, scale*lat)
        plt.plot(x, y, 'g')

    for longitude in range(-limit, limit+1, step):
        x, y = guyou(scale*longitude, scale*long_line)
        plt.plot(x, y, 'b')
    #plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.title('forward transform')

    plt.subplot(212)
    for lat in range(-limit, limit+1, step):
        x, y = guyou_invert(scale*lat_line, scale*lat)
        plt.plot(x, y, 'g')

    for long in range(-limit, limit+1, step):
        x, y = guyou_invert(scale*long, scale*long_line)
        plt.plot(x, y, 'b')
    plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.title('inverse transform')

def main():
    """Show the Guyou transformation"""
    plot_grid()
    import matplotlib.pyplot as plt
    plt.show()

if __name__ == "__main__":
    main()

_ = """
// Javascript source for elliptic functions
//
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

  // [PAK] Cephes uses dn = sqrt(1 - m*sin^2 phi) rather than cos(phi)/cos(phi-b)
  // DLMF says the second version is unstable near x = K.
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


export function guyouRaw(lambda, phi) {
  var k_ = (sqrt2 - 1) / (sqrt2 + 1),
      k = sqrt(1 - k_ * k_),
      K = ellipticF(halfPi, k * k),
      f = -1,
      psi = log(tan(pi / 4 + abs(phi) / 2)),
      r = exp(f * psi) / sqrt(k_),
      at = guyouComplexAtan(r * cos(f * lambda), r * sin(f * lambda)),
      t = ellipticFi(at[0], at[1], k * k);
  return [-t[1], (phi >= 0 ? 1 : -1) * (0.5 * K - t[0])];
}

function guyouComplexAtan(x, y) {
  var x2 = x * x,
      y_1 = y + 1,
      t = 1 - x2 - y * y;
  return [
   0.5 * ((x >= 0 ? halfPi : -halfPi) - atan2(t, 2 * x)),
    -0.25 * log(t * t + 4 * x2) +0.5 * log(y_1 * y_1 + x2)
  ];
}

function guyouComplexDivide(a, b) {
  var denominator = b[0] * b[0] + b[1] * b[1];
  return [
    (a[0] * b[0] + a[1] * b[1]) / denominator,
    (a[1] * b[0] - a[0] * b[1]) / denominator
  ];
}

guyouRaw.invert = function(x, y) {
  var k_ = (sqrt2 - 1) / (sqrt2 + 1),
      k = sqrt(1 - k_ * k_),
      K = ellipticF(halfPi, k * k),
      f = -1,
      j = ellipticJi(0.5 * K - y, -x, k * k),
      tn = guyouComplexDivide(j[0], j[1]),
      lambda = atan2(tn[1], tn[0]) / f;
  return [
    lambda,
    2 * atan(exp(0.5 / f * log(k_ * tn[0] * tn[0] + k_ * tn[1] * tn[1]))) - halfPi
  ];
};
"""
