double Iq(double q,
    double phiA, double nA, double vA, double lA, double bA,
    double nB, double vB, double lB, double bB,
    double chiAB
    );

static inline double debye(double X) {
  if (X < 1e-3) {
    // 2*(e^-X - 1 + X)/X^2 ≈ 1 - X/3 + X^2/12 - X^3/60 + ...
    const double X2 = X*X;
    return 1.0 - X/3.0 + X2/12.0 - X2*X/60.0;
  } else {
    // use expm1 for better accuracy when X is small-to-moderate
    return 2.0*(expm1(-X) + X)/ (X*X);
  }
}

double Iq(double q,
    double phiA,  // VOL FRACTION
    double nA,    // DEGREE OF POLYMERIZATION
    double vA,    // SPECIFIC VOLUME
    double lA,    // SCATT. LENGTH
    double bA,    // SEGMENT LENGTH
    double nB,    // DEGREE OF POLYMERIZATION
    double vB,    // SPECIFIC VOLUME
    double lB,    // SCATT. LENGTH
    double bB,    // SEGMENT LENGTH
    double chiAB  // CHI PARAM
    )
{
  const double q2 = q*q;

  // Volume fraction of polymer B
  const double phiB = 1.0 - phiA;

  // Calculate Q^2 * Rg^2 for each homopolymer assuming random walk
  const double XA = q2 * bA * bA * nA / 6.0;
  const double XB = q2 * bB * bB * nB / 6.0;

  //calculate all partial structure factors Pij and normalize n^2
  const double PAA = debye(XA);
  const double PBB = debye(XB);
  const double S0AA = nA * phiA * vA * PAA;
  const double S0BB = nB * phiB * vB *PBB;

  const double v0 = sqrt(vA*vB); // Reference volume

  const double v11 = (1.0/S0BB) - (2.0*chiAB/v0);

  const double Nav = 6.022141e+23;
  const double contrast = (lA/vA -lB/vB) * (lA/vA -lB/vB) * Nav;

  double intensity = contrast * S0AA / (1.0 + v11*S0AA);

  //rescale for units of lA^2 and lB^2 (fm^2 to cm^2)
  intensity *= 1.0e-26;

  return intensity;
}
