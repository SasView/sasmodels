// acos(-1.0/3.0)/2.0 = half of tetrahedral angle ~54.7356 deg
#define TETRAHEDRAL_HALF_ANGLE acos(-1.0 / 3.0) / 2.0

static double u_n(int n, double theta, double alpha) {
  const double phi[4] = {0.0, M_PI_2, M_PI, 3.0 * M_PI_2};
  const double sign[4] = {1.0, -1.0, 1.0, -1.0};
  return sign[n] * cos(TETRAHEDRAL_HALF_ANGLE) * cos(theta) +
         sin(TETRAHEDRAL_HALF_ANGLE) * sin(theta) * cos(alpha - phi[n]);
}

// L: arm length, R: outer radius, t: shell thickness, R_c = R - t: core radius
static double Fq_n(double q, double u, double L, double R, double t,
                   double contrast_core, double contrast_shell) {
  double R_c = R - t;
  double quL2 = q * u * L * 0.5;
  double mu = sqrt(fmax(0.0, 1.0 - u * u));
  double V_c = M_PI * R_c * R_c * L;
  double V = M_PI * R * R * L;
  return sas_sinx_x(quL2) *
         (contrast_core * V_c * sas_2J1x_x(q * mu * R_c) +
          contrast_shell * V * sas_2J1x_x(q * mu * R));
}

static double form_volume(double L, double R, double t) {
  // V = 4 * pi * R^2 * L  (outer volume of 4 arms)
  return 4.0 * M_PI * R * R * L;
}

static double Iq(double q, double L, double R, double t,
                 double sld_core, double sld_shell, double sld_solvent) {
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

      double sum_arms = 0.0;
      for (int n = 0; n < 4; n++) {
        for (int m = 0; m < 4; m++) {
          double u = u_n(n, theta, alpha);
          double v = u_n(m, theta, alpha);
          double Fn = Fq_n(q, u, L, R, t, contrast_core, contrast_shell);
          double Fm = Fq_n(q, v, L, R, t, contrast_core, contrast_shell);
          sum_arms += Fn * Fm * cos(q * (u - v) * L / 2.0) * sin(theta);
        }
      }
      integral_alpha += sum_arms * w_alpha;
    }
    total += integral_alpha * w_theta;
  }
  return 1e-4 * total / (4.0 * M_PI);
}
