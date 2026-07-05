// Half of the tetrahedral angle acos(-1/3)/2 ~54.7356 deg.
// Precompute its cosine and sine: cos = sqrt(1/3), sin = sqrt(2/3).
#define COS_HALF_ANGLE 0.57735026918962573  // sqrt(1/3)
#define SIN_HALF_ANGLE 0.81649658092772603  // sqrt(2/3)

static double u_n(int n, double theta, double alpha) {
  const double phi[4] = {0.0, M_PI_2, M_PI, 3.0 * M_PI_2};
  const double sign[4] = {1.0, -1.0, 1.0, -1.0};
  return sign[n] * COS_HALF_ANGLE * cos(theta) +
         SIN_HALF_ANGLE * sin(theta) * cos(alpha - phi[n]);
}

// L: arm length, R: core radius, t: shell thickness, R_o = R + t: outer radius
static double Fq_n(double q, double u, double L, double R, double t,
                   double contrast_core, double contrast_shell) {
  double R_o = R + t;
  double quL2 = q * u * L * 0.5;
  double mu = sqrt(fmax(0.0, 1.0 - u * u));
  double V_c = M_PI * R * R * L;
  double V_o = M_PI * R_o * R_o * L;
  return sas_sinx_x(quL2) * (contrast_core * V_c * sas_2J1x_x(q * mu * R) +
                             contrast_shell * V_o * sas_2J1x_x(q * mu * R_o));
}

static double form_volume(double L, double R, double t) {
  // V = 4 * pi * (R + t)^2 * L  (outer volume of 4 arms)
  return 4.0 * M_PI * (R + t) * (R + t) * L;
}

static double radius_effective(int mode, double L, double R, double t) {
  switch (mode) {
    default:
    case 1:  // equivalent volume sphere
      return cbrt(form_volume(L, R, t) / M_4PI_3);
    case 2:  // length of tetrapod arms (L)
      return L;
  }
}

static double Iq(double q, double L, double R, double t, double sld_core,
                 double sld_shell, double sld_solvent) {
  double contrast_core = sld_core - sld_shell;
  double contrast_shell = sld_shell - sld_solvent;
  double total = 0.0;

  for (int dtheta = 0; dtheta < GAUSS_N; dtheta++) {
    double theta =
        M_PI_2 * (GAUSS_Z[dtheta] + 1.0);  // map from [-1, 1] to [0, pi]
    double w_theta =
        GAUSS_W[dtheta] * M_PI_2;  // adjust weight for the new range

    double integral_alpha = 0.0;
    for (int dalpha = 0; dalpha < GAUSS_N; dalpha++) {
      double alpha =
          M_PI * (GAUSS_Z[dalpha] + 1.0);  // map from [-1, 1] to [0, 2*pi]
      double w_alpha =
          GAUSS_W[dalpha] * M_PI;  // adjust weight for the new range

      double u[4], F[4];
      for (int n = 0; n < 4; n++) {
        u[n] = u_n(n, theta, alpha);
        F[n] = Fq_n(q, u[n], L, R, t, contrast_core, contrast_shell);
      }
      double sum_arms = 0.0;
      for (int n = 0; n < 4; n++) {
        sum_arms += F[n] * F[n];
        for (int m = n + 1; m < 4; m++) {
          sum_arms += 2.0 * F[n] * F[m] * cos(q * (u[n] - u[m]) * L / 2.0);
        }
      }
      sum_arms *= sin(theta);
      integral_alpha += sum_arms * w_alpha;
    }
    total += integral_alpha * w_theta;
  }
  return 1e-4 * total / (4.0 * M_PI);
}
