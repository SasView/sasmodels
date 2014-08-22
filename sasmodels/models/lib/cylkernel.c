// Compute the form factor for a cylinder
// qq is the q-value for the calculation (1/A)
// radius is the radius of the cylinder (A)
// h is the HALF-LENGTH of the cylinder = L/2 (A)
real CylKernel(real q, real radius, real h, real theta);
real CylKernel(real q, real radius, real h, real theta)
{
    real sn, cn;
    SINCOS(theta, sn, cn);
    const real besarg = q*radius*sn;
    const real siarg = q*h*cn;
    // lim_{x->0} J1(x)/x = 1/2,  lim_{x->0} sin(x)/x = 1
    const real bj = (besarg == REAL(0.0) ? REAL(0.5) : J1(besarg)/besarg);
    const real si = (siarg == REAL(0.0) ? REAL(1.0) : sin(siarg)/siarg);
    return REAL(4.0)*sin(theta)*bj*bj*si*si;
}
