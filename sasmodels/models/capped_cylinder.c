real form_volume(real radius, real cap_radius, real length);
real Iq(real q, real sld, real solvent_sld,
    real radius, real cap_radius, real length);
real Iqxy(real qx, real qy, real sld, real solvent_sld,
    real radius, real cap_radius, real length, real theta, real phi);

// Integral over a convex lens kernel for t in [h/R,1].  See the docs for
// the definition of the function being integrated.
//   q is the magnitude of the q vector.
//   h is the length of the lens "inside" the cylinder.  This negative wrt the
//       definition of h in the docs.
//   cap_radius is the radius of the lens
//   length is the cylinder length, or the separation between the lens halves
//   alpha is the angle of the cylinder wrt q.
real _cap_kernel(real q, real h, real cap_radius, real length,
                 real sin_alpha, real cos_alpha);
real _cap_kernel(real q, real h, real cap_radius, real length,
                 real sin_alpha, real cos_alpha)
{
    // For speed, we are pre-calculating terms which are constant over the
    // kernel.
    const real upper = REAL(1.0);
    const real lower = h/cap_radius; // integral lower bound
    // cos term in integral is:
    //    cos (q (R t - h + L/2) cos(alpha))
    // so turn it into:
    //    cos (m t + b)
    // where:
    //    m = q R cos(alpha)
    //    b = q(L/2-h) cos(alpha)
    const real m = q*cap_radius*cos_alpha; // cos argument slope
    const real b = q*(REAL(0.5)*length-h)*cos_alpha; // cos argument intercept
    const real qrst = q*cap_radius*sin_alpha; // Q*R*sin(theta)
    real total = REAL(0.0);
    for (int i=0; i<76 ;i++) {
        // translate a point in [-1,1] to a point in [lower,upper]
        //const real t = ( Gauss76Z[i]*(upper-lower) + upper + lower )/2.0;
        const real t = REAL(0.5)*(Gauss76Z[i]*(upper-lower)+upper+lower);
        const real radical = REAL(1.0) - t*t;
        const real arg = qrst*sqrt(radical); // cap bessel function arg
        const real be = (arg == REAL(0.0) ? REAL(0.5) : J1(arg)/arg);
        const real Fq = cos(m*t + b) * radical * be;
        total += Gauss76Wt[i] * Fq;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    //const real form = (upper-lower)/2.0*total;
    const real integral = REAL(0.5)*(upper-lower)*total;
    return REAL(4.0)*M_PI*cap_radius*cap_radius*cap_radius*integral;
}

real form_volume(real radius, real cap_radius, real length)
{
    // cap radius should never be less than radius when this is called
    // Note: cap volume = pi hc/6 * (3 a^2 + hc^2), where a is the cylinder
    // radius and hc is the height of the cap.  Multiply by two for both ends.
    // So:
    //      cap V = pi hc (r^2 + hc^2/3)
    //      cylinder V = pi r^2 L
    //      V = cylinder V + cap V
    //        = pi r^2 L + pi hc (r^2 + hc^2/3)
    //        = pi * (r^2 (L+hc) + hc^3/3)
    const real hc = cap_radius - sqrt(cap_radius*cap_radius - radius*radius);
    return M_PI*(radius*radius*(length+hc) + REAL(0.333333333333333)*hc*hc*hc);
}

real Iq(real q,
    real sld,
    real solvent_sld,
    real radius,
    real cap_radius,
    real length)
{
    real sn, cn; // slots to hold sincos function output

    // Exclude invalid inputs.
    if (cap_radius < radius) return REAL(-1.0);

    const real lower = REAL(0.0);
    const real upper = M_PI_2;
    const real h = sqrt(cap_radius*cap_radius - radius*radius); // negative h
    real total = REAL(0.0);
    for (int i=0; i<76 ;i++) {
        // translate a point in [-1,1] to a point in [lower,upper]
        const real alpha= REAL(0.5)*(Gauss76Z[i]*(upper-lower) + upper + lower);
        SINCOS(alpha, sn, cn);

        const real cap_Fq = _cap_kernel(q, h, cap_radius, length, sn, cn);

        // The following is CylKernel() / sin(alpha), but we are doing it in place
        // to avoid sin(alpha)/sin(alpha) for alpha = 0.  It is also a teensy bit
        // faster since we don't multiply and divide sin(alpha).
        const real besarg = q*radius*sn;
        const real siarg = q*REAL(0.5)*length*cn;
        // lim_{x->0} J1(x)/x = 1/2,   lim_{x->0} sin(x)/x = 1
        const real bj = (besarg == REAL(0.0) ? REAL(0.5) : J1(besarg)/besarg);
        const real si = (siarg == REAL(0.0) ? REAL(1.0) : sin(siarg)/siarg);
        const real cyl_Fq = M_PI*radius*radius*length*REAL(2.0)*bj*si;

        // Volume weighted average F(q)
        const real Aq = cyl_Fq + cap_Fq;
        total += Gauss76Wt[i] * Aq * Aq * sn; // sn for spherical coord integration
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    const real form = total * REAL(0.5)*(upper-lower);

    // Multiply by contrast^2, normalize by cylinder volume and convert to cm-1
    // NOTE that for this (Fournet) definition of the integral, one must MULTIPLY by Vcyl
    // The additional volume factor is for polydisperse volume normalization.
    const real s = (sld - solvent_sld);
    return REAL(1.0e-4) * form * s * s; // form_volume(radius, cap_radius, length);
}


real Iqxy(real qx, real qy,
    real sld,
    real solvent_sld,
    real radius,
    real cap_radius,
    real length,
    real theta,
    real phi)
{
    real sn, cn; // slots to hold sincos function output

    // Exclude invalid inputs.
    if (cap_radius < radius) return REAL(-1.0);

    // Compute angle alpha between q and the cylinder axis
    SINCOS(theta*M_PI_180, sn, cn);
    // # The following correction factor exists in sasview, but it can't be
    // # right, so we are leaving it out for now.
    const real q = sqrt(qx*qx+qy*qy);
    const real cos_val = cn*cos(phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const real alpha = acos(cos_val); // rod angle relative to q
    SINCOS(alpha, sn, cn);

    const real h = sqrt(cap_radius*cap_radius - radius*radius); // negative h
    const real cap_Fq = _cap_kernel(q, h, cap_radius, length, sn, cn);

    const real besarg = q*radius*sn;
    const real siarg = q*REAL(0.5)*length*cn;
    // lim_{x->0} J1(x)/x = 1/2,   lim_{x->0} sin(x)/x = 1
    const real bj = (besarg == REAL(0.0) ? REAL(0.5) : J1(besarg)/besarg);
    const real si = (siarg == REAL(0.0) ? REAL(1.0) : sin(siarg)/siarg);
    const real cyl_Fq = M_PI*radius*radius*length*REAL(2.0)*bj*si;

    // Volume weighted average F(q)
    const real Aq = cyl_Fq + cap_Fq;

    // Multiply by contrast^2, normalize by cylinder volume and convert to cm-1
    const real s = (sld - solvent_sld);
    return REAL(1.0e-4) * Aq * Aq * s * s; // form_volume(radius, cap_radius, length);
}
