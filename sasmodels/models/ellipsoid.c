real form_volume(real rpolar, real requatorial);
real Iq(real q, real sld, real solvent_sld, real rpolar, real requatorial);
real Iqxy(real qx, real qy, real sld, real solvent_sld,
    real rpolar, real requatorial, real theta, real phi);

real _ellipsoid_kernel(real q, real rpolar, real requatorial, real cos_alpha);
real _ellipsoid_kernel(real q, real rpolar, real requatorial, real cos_alpha)
{
    real sn, cn;
    real ratio = rpolar/requatorial;
    const real u = q*requatorial*sqrt(REAL(1.0)
                   + cos_alpha*cos_alpha*(ratio*ratio - REAL(1.0)));
    SINCOS(u, sn, cn);
    const real f = ( u==REAL(0.0) ? REAL(1.0) : REAL(3.0)*(sn-u*cn)/(u*u*u) );
    return f*f;
}

real form_volume(real rpolar, real requatorial)
{
    return REAL(1.333333333333333)*M_PI*rpolar*requatorial*requatorial;
}

real Iq(real q,
    real sld,
    real solvent_sld,
    real rpolar,
    real requatorial)
{
    //const real lower = REAL(0.0);
    //const real upper = REAL(1.0);
    real total = REAL(0.0);
    for (int i=0;i<76;i++) {
        //const real cos_alpha = (Gauss76Z[i]*(upper-lower) + upper + lower)/2;
        const real cos_alpha = REAL(0.5)*(Gauss76Z[i] + REAL(1.0));
        total += Gauss76Wt[i] * _ellipsoid_kernel(q, rpolar, requatorial, cos_alpha);
    }
    //const real form = (upper-lower)/2*total;
    const real form = REAL(0.5)*total;
    const real s = (sld - solvent_sld) * form_volume(rpolar, requatorial);
    return REAL(1.0e-4) * form * s * s;
}

real Iqxy(real qx, real qy,
    real sld,
    real solvent_sld,
    real rpolar,
    real requatorial,
    real theta,
    real phi)
{
    real sn, cn;

    const real q = sqrt(qx*qx + qy*qy);
    SINCOS(theta*M_PI_180, sn, cn);
    const real cos_alpha = cn*cos(phi*M_PI_180)*(qx/q) + sn*(qy/q);
    const real form = _ellipsoid_kernel(q, rpolar, requatorial, cos_alpha);
    const real s = (sld - solvent_sld) * form_volume(rpolar, requatorial);

    return REAL(1.0e-4) * form * s * s;
}

