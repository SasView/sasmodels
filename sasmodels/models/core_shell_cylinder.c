real form_volume(real radius, real thickness, real length);
real Iq(real q, real core_sld, real shell_sld, real solvent_sld,
    real radius, real thickness, real length);
real Iqxy(real qx, real qy, real core_sld, real shell_sld, real solvent_sld,
    real radius, real thickness, real length, real theta, real phi);

// twovd = 2 * volume * delta_rho
// besarg = q * R * sin(alpha)
// siarg = q * L/2 * cos(alpha)
real _cyl(real twovd, real besarg, real siarg);
real _cyl(real twovd, real besarg, real siarg)
{
    const real bj = (besarg == REAL(0.0) ? REAL(0.5) : J1(besarg)/besarg);
    const real si = (siarg == REAL(0.0) ? REAL(1.0) : sin(siarg)/siarg);
    return twovd*si*bj;
}

real form_volume(real radius, real thickness, real length)
{
    return M_PI*(radius+thickness)*(radius+thickness)*(length+2*thickness);
}

real Iq(real q,
    real core_sld,
    real shell_sld,
    real solvent_sld,
    real radius,
    real thickness,
    real length)
{
    // precalculate constants
    const real core_qr = q*radius;
    const real core_qh = q*REAL(0.5)*length;
    const real core_twovd = REAL(2.0) * form_volume(radius,0,length)
                            * (core_sld-shell_sld);
    const real shell_qr = q*(radius + thickness);
    const real shell_qh = q*(REAL(0.5)*length + thickness);
    const real shell_twovd = REAL(2.0) * form_volume(radius,thickness,length)
                             * (shell_sld-solvent_sld);
    real total = REAL(0.0);
    // real lower=0, upper=M_PI_2;
    for (int i=0; i<76 ;i++) {
        // translate a point in [-1,1] to a point in [lower,upper]
        //const real alpha = ( Gauss76Z[i]*(upper-lower) + upper + lower )/2.0;
        real sn, cn;
        const real alpha = REAL(0.5)*(Gauss76Z[i]*M_PI_2 + M_PI_2);
        SINCOS(alpha, sn, cn);
        const real fq = _cyl(core_twovd, core_qr*sn, core_qh*cn)
            + _cyl(shell_twovd, shell_qr*sn, shell_qh*cn);
        total += Gauss76Wt[i] * fq * fq * sn;
    }
    // translate dx in [-1,1] to dx in [lower,upper]
    //const real form = (upper-lower)/2.0*total;
    return REAL(1.0e-4) * total * M_PI_4;
}


real Iqxy(real qx, real qy,
    real core_sld,
    real shell_sld,
    real solvent_sld,
    real radius,
    real thickness,
    real length,
    real theta,
    real phi)
{
    real sn, cn; // slots to hold sincos function output

    // Compute angle alpha between q and the cylinder axis
    SINCOS(theta*M_PI_180, sn, cn);
    // # The following correction factor exists in sasview, but it can't be
    // # right, so we are leaving it out for now.
    // const real correction = fabs(cn)*M_PI_2;
    const real q = sqrt(qx*qx+qy*qy);
    const real cos_val = cn*cos(phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const real alpha = acos(cos_val);

    const real core_qr = q*radius;
    const real core_qh = q*REAL(0.5)*length;
    const real core_twovd = REAL(2.0) * form_volume(radius,0,length)
                            * (core_sld-shell_sld);
    const real shell_qr = q*(radius + thickness);
    const real shell_qh = q*(REAL(0.5)*length + thickness);
    const real shell_twovd = REAL(2.0) * form_volume(radius,thickness,length)
                             * (shell_sld-solvent_sld);

    SINCOS(alpha, sn, cn);
    const real fq = _cyl(core_twovd, core_qr*sn, core_qh*cn)
        + _cyl(shell_twovd, shell_qr*sn, shell_qh*cn);
    return REAL(1.0e-4) * fq * fq;
}
