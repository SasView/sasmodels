
static double
f_exp(double q, double r, double sld_in, double sld_out,
    double thickness, double A)
{
  const double vol = M_4PI_3 * cube(r);
  const double qr = q * r;
  const double alpha = A * r/thickness;
  const double bes = sph_j1c(qr);
  const double B = (sld_out - sld_in)/expm1(A);
  const double C = sld_in - B;
  double fun;
  if (qr == 0.0) {
    fun = 1.0;
  } else if (fabs(A) > 0.0) {
    const double qrsq = qr*qr;
    const double alphasq = alpha*alpha;
    const double sumsq = alphasq + qrsq;
    double sinqr, cosqr;
    SINCOS(qr, sinqr, cosqr);
    fun = -3.0*(
            ((alphasq - qrsq)*sinqr/qr - 2.0*alpha*cosqr) / sumsq
                - (alpha*sinqr/qr - cosqr)
        ) / sumsq;
  } else {
    fun = bes;
  }
  return vol * (B*fun + C*bes);
}

static double
f_linear(double q, double r, double sld, double slope)
{
  const double vol = M_4PI_3 * cube(r);
  const double qr = q * r;
  const double bes = sph_j1c(qr);
  double fun = 0.0;
  if (qr > 0.0) {
    const double qrsq = qr*qr;
    double sinqr, cosqr;
    SINCOS(qr, sinqr, cosqr);
    // Jae-He's code seems to simplify to this
    //     fun = 3.0 * slope * r * (2.0*qr*sinqr - (qrsq-2.0)*cosqr)/(qrsq*qrsq);
    // Rederiving the math, we get the following instead:
    fun = 3.0 * slope * r * (2.0*cosqr + qr*sinqr)/(qrsq*qrsq);
  }
  return vol * (sld*bes + fun);
}

static double
f_constant(double q, double r, double sld)
{
  const double bes = sph_j1c(q * r);
  const double vol = M_4PI_3 * cube(r);
  return sld * vol * bes;
}

static double
form_volume(double core_radius, double n, double thickness[])
{
  int i;
  double r = core_radius;
  for (i=0; i < n; i++) {
    r += thickness[i];
  }
  return M_4PI_3*cube(r);
}

static double
Iq(double q, double sld_core, double core_radius, double sld_solvent,
    double n, double sld_in[], double sld_out[], double thickness[],
    double A[])
{
  int i;
  double r = core_radius;
  double f = f_constant(q, r, sld_core);
  for (i=0; i<n; i++){
    const double r0 = r;
    r += thickness[i];
    if (r == r0) {
      // no thickness, so nothing to add
    } else if (fabs(A[i]) < 1e-16 || sld_out[i] == sld_in[i]) {
      f -= f_constant(q, r0, sld_in[i]);
      f += f_constant(q, r, sld_in[i]);
    } else if (fabs(A[i]) < 1e-4) {
      const double slope = (sld_out[i] - sld_in[i])/thickness[i];
      f -= f_linear(q, r0, sld_in[i], slope);
      f += f_linear(q, r, sld_out[i], slope);
    } else {
      f -= f_exp(q, r0, sld_in[i], sld_out[i], thickness[i], A[i]);
      f += f_exp(q, r, sld_in[i], sld_out[i], thickness[i], A[i]);
    }
  }
  f -= f_constant(q, r, sld_solvent);
  const double f2 = f * f * 1.0e-4;

  return f2;
}
