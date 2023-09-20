// from https://raw.githubusercontent.com/opencv/opencv/f6c573880e57e78e0c07d8832ee6015f2337a019/modules/calib3d/src/polynom_solver.cpp
// with serious modifications July 2019 for sasview
/*
By downloading, copying, installing or using the software you agree to this license.
If you do not agree to this license, do not download, install,
copy or use the software.


                          License Agreement
               For Open Source Computer Vision Library
                       (3-clause BSD License)

Copyright (C) 2000-2019, Intel Corporation, all rights reserved.
Copyright (C) 2009-2011, Willow Garage Inc., all rights reserved.
Copyright (C) 2009-2016, NVIDIA Corporation, all rights reserved.
Copyright (C) 2010-2013, Advanced Micro Devices, Inc., all rights reserved.
Copyright (C) 2015-2016, OpenCV Foundation, all rights reserved.
Copyright (C) 2015-2016, Itseez Inc., all rights reserved.
Third party copyrights are property of their respective owners.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * Neither the names of the copyright holders nor the names of the contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

This software is provided by the copyright holders and contributors "as is" and
any express or implied warranties, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose are disclaimed.
In no event shall copyright holders or contributors be liable for any direct,
indirect, incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused
and on any theory of liability, whether in contract, strict liability,
or tort (including negligence or otherwise) arising in any way out of
the use of this software, even if advised of the possibility of such damage.
*/

int solve_deg2(double a, double b, double c, double *x1, double *x2);
int solve_deg2(double a, double b, double c, double *x1, double *x2)
{
  double delta = b * b - 4 * a * c;

  if (delta < 0) return 0;

  double inv_2a = 0.5 / a;

  if (delta == 0) {
    *x1 = *x2 = -b * inv_2a;
    return 1;
  }

  double sqrt_delta = sqrt(delta);
  *x1 = (-b + sqrt_delta) * inv_2a;
  *x2 = (-b - sqrt_delta) * inv_2a;
  return 2;
}


/// Reference : Eric W. Weisstein. "Cubic Equation." From MathWorld--A Wolfram Web Resource.
/// http://mathworld.wolfram.com/CubicEquation.html
///  Modified to return thick_shell, the largest root of bespoke equation

double
solve_shell_v3(double radius_equat_core,
    double x_core,
    double vol_dry_shell_over_core,
    double x_polar_shell,
    double f_solvent_in_shell)
{
    double x0,x1,x2;
    int nroots;
    const double a = x_polar_shell;
    const double b = radius_equat_core*(x_core + 2.0*x_polar_shell);
    const double c = square(radius_equat_core)*(2.0*x_core + x_polar_shell);
    const double d = cube(radius_equat_core)*x_core*vol_dry_shell_over_core/(f_solvent_in_shell -1.0);
  if (a == 0) {
    // Solve second order system
    if (b == 0)	{
      // Solve first order system
      if (c == 0)
        return 0.0;
      x0 = -d / c;
      //printf("quadratic  x0 = %g\n",x0);
      return x0;
    }
    nroots = solve_deg2(b, c, d, &x0, &x1);
    if(nroots>1){
          if(x1 > x0){x0 = x1;}}
          //printf("quadratic x0 = %g\n",x0);
    return x0;
  }
  // Calculate the normalized form x^3 + a2 * x^2 + a1 * x + a0 = 0
  double inv_a = 1. / a;
  double b_a = inv_a * b, b_a2 = b_a * b_a;
  double c_a = inv_a * c;
  double d_a = inv_a * d;

  // Solve the cubic equation
  double Q = (3 * c_a - b_a2) / 9;
  double R = (9 * b_a * c_a - 27 * d_a - 2 * b_a * b_a2) / 54;
  double Q3 = Q * Q * Q;
  double D = Q3 + R * R;
  double b_a_3 = (1. / 3.) * b_a;

  if (Q == 0) {
    if(R == 0) {
          // three identical roots
      x0 = - b_a_3;
    } else {
      x0 = copysign(cbrt(fabs(2.0*R)), R) - b_a_3;
    } }
  else {
  if (D <= 0) {
    // Three real roots
    double theta = acos(R / sqrt(-Q3));
    double sqrt_Q = sqrt(-Q);
    x0 = 2 * sqrt_Q * cos(theta             / 3.0) - b_a_3;
    x1 = 2 * sqrt_Q * cos((theta + 2 * M_PI)/ 3.0) - b_a_3;
    x2 = 2 * sqrt_Q * cos((theta + 4 * M_PI)/ 3.0) - b_a_3;
    if(x2 > x1){x1 = x2;}
    if(x1 > x0){x0 = x1;}
  }  else {
  // D > 0, only one real root
  //double AD = pow(fabs(R) + sqrt(D), 1.0 / 3.0) * (R > 0 ? 1 : (R < 0 ? -1 : 0));
  double AD =  copysign(cbrt(fabs(R) + sqrt(D)), R+sqrt(D));
  double BD = (AD == 0) ? 0 : -Q / AD;
  // Calculate the only real root
  x0 = AD + BD - b_a_3;
  } }
  // improve solution with two Newton-Raphson iterations, as rounding errors often occur in the algebra
  // (In most cases five or six of these would find the thick_shell solution given a reasonable starting guess ...)
  x0 = x0 - (((x0 + b_a)*x0 + c_a)*x0 + d_a)/((3.0*x0 + 2.0*b_a)*x0 + c_a);
  x0 = x0 - (((x0 + b_a)*x0 + c_a)*x0 + d_a)/((3.0*x0 + 2.0*b_a)*x0 + c_a);
  //printf(" x0 = %g\n",x0);
  return x0;
}


