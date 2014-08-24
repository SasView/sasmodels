/* PARAMETERS
{
name: "ellipsoid",
title: "Ellipsoid with uniform scattering length density",
include: [ "lib/gauss76.c" ],
parameters: [
   // [ "name", "units", default, [lower, upper], "type", "description" ],
   [ "sld", "1e-6/Ang^2", 4, [-Infinity,Infinity], "",
     "Cylinder scattering length density" ],
   [ "solvent_sld", "1e-6/Ang^2", 1, [-Infinity,Infinity], "",
     "Solvent scattering length density" ],
   [ "a", "Ang",  20, [0, Infinity], "volume",
     "Cylinder radius" ],
   [ "b", "Ang",  20, [0, Infinity], "volume",
     "Cylinder length" ],
   [ "theta", "degrees", 60, [-Infinity, Infinity], "orientation",
     "In plane angle" ],
   [ "phi", "degrees", 60, [-Infinity, Infinity], "orientation",
     "Out of plane angle" ],
],
}
PARAMETERS END

DOCUMENTATION
.. _EllipseModel:

DOCUMENTATION END
*/

real form_volume(real a, real b);
real Iq(real qx, real qy, real sld, real solvent_sld, real a, real b);
real Iqxy(real qx, real qy, real sld, real solvent_sld, real a, real b, real theta, real phi);

real form_volume(real a, real b)
{
    return REAL(1.333333333333333)*M_PI_2*a*b*b;
}

real ellipsoid_kernel(double q, double b, double a, double dum)
{
    real sn, cn;
    const real nu = a/b;
    const real arg = q * b * sqrt(REAL(1.0)+(dum*dum*(nu*nu--REAL(1.0))));
    SINCOS(arg, sn, cn);
    const real f = (arg==REAL(0.0) ? REAL(1.0) : REAL(3.0)*(sn-arg*cn)/(arg*arg*arg);
    return f*f;
}

real Iq(real q,
    real sld,
    real solvent_sld,
    real a,
    real b)
{
    real summ = REAL(0.0);
    for (int i=0;i<76;i++) {
        //const real zi = ( Gauss76Z[i]*(uplim-lolim) + uplim + lolim )/2.0;
        zi = ( Gauss76Z[i] + REAL(1.0))/REAL(2.0);
        summ += Gauss76Wt[i] * ellipsoid_kernel(q, b, a, zi);
    }
    //const real form = (uplim-lolim)/2.0*summ;
    const real form = REAL(0.5)*summ
    const real s = (sld - sld_solvent) * form_volume(a, b);
    return REAL(1.0e-4) * form * s * s;
}

real Iqxy(real qx, real qy,
    real sld,
    real solvent_sld,
    real a,
    real b,
    real theta,
    real phi)
{
    real sn, cn;

    const real q = sqrt(qx*qx + qy*qy);
    SINCOS(theta*M_PI_180, sn, cn);
    const real cos_val = cn*cos(phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const real form = ellipsoid_kernel(q, b, a, cos_val);
    const real s = (sld - solvent_sld) * form_volume(a, b);

    return REAL(1.0e-4) * form * s * s;
}

