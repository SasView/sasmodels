real form_volume(real req_minor, real req_major, real rpolar);
real Iq(real q, real sld, real solvent_sld,
    real req_minor, real req_major, real rpolar);
real Iqxy(real qx, real qy, real sld, real solvent_sld,
    real req_minor, real req_major, real rpolar, real theta, real phi, real psi);

real form_volume(real req_minor, real req_major, real rpolar)
{
    return REAL(1.333333333333333)*M_PI*req_minor*req_major*rpolar;
}

real Iq(real q,
    real sld,
    real solvent_sld,
    real req_minor,
    real req_major,
    real rpolar)
{
    // if (req_minor > req_major || req_major > rpolar) return REAL(-1.0);  // Exclude invalid region

    real sn, cn;
    real st, ct;
    //const real lower = REAL(0.0);
    //const real upper = REAL(1.0);
    real outer = REAL(0.0);
    for (int i=0;i<76;i++) {
        //const real cos_alpha = (Gauss76Z[i]*(upper-lower) + upper + lower)/2;
        const real x = REAL(0.5)*(Gauss76Z[i] + REAL(1.0));
        SINCOS(M_PI_2*x, sn, cn);
        const real acosx2 = req_minor*req_minor*cn*cn;
        const real bsinx2 = req_major*req_major*sn*sn;
        const real c2 = rpolar*rpolar;

        real inner = REAL(0.0);
        for (int j=0;j<76;j++) {
            const real y = REAL(0.5)*(Gauss76Z[j] + REAL(1.0));
            const real t = q*sqrt(acosx2 + bsinx2*(REAL(1.0)-y*y) + c2*y*y);
            SINCOS(t, st, ct);
            const real fq = ( t==REAL(0.0) ? REAL(1.0) : REAL(3.0)*(st-t*ct)/(t*t*t) );
            inner += Gauss76Wt[j] * fq * fq ;
        }
        outer += Gauss76Wt[i] * REAL(0.5) * inner;
    }
    //const real fq2 = (upper-lower)/2*outer;
    const real fq2 = REAL(0.5)*outer;
    const real s = (sld - solvent_sld) * form_volume(req_minor, req_major, rpolar);
    return REAL(1.0e-4) * fq2 * s * s;
}

real Iqxy(real qx, real qy,
    real sld,
    real solvent_sld,
    real req_minor,
    real req_major,
    real rpolar,
    real theta,
    real phi,
    real psi)
{
    // if (req_minor > req_major || req_major > rpolar) return REAL(-1.0);  // Exclude invalid region

    real stheta, ctheta;
    real sphi, cphi;
    real spsi, cpsi;
    real st, ct;

    const real q = sqrt(qx*qx + qy*qy);
    const real qxhat = qx/q;
    const real qyhat = qy/q;
    SINCOS(theta*M_PI_180, stheta, ctheta);
    SINCOS(phi*M_PI_180, sphi, cphi);
    SINCOS(psi*M_PI_180, spsi, cpsi);
    const real calpha = ctheta*cphi*qxhat + stheta*qyhat;
    const real cnu = (-cphi*spsi*stheta + sphi*cpsi)*qxhat + spsi*ctheta*qyhat;
    const real cmu = (-stheta*cpsi*cphi - spsi*sphi)*qxhat + ctheta*cpsi*qyhat;
    const real t = q*sqrt(req_minor*req_minor*cnu*cnu
                          + req_major*req_major*cmu*cmu
                          + rpolar*rpolar*calpha*calpha);
    SINCOS(t, st, ct);
    const real fq = ( t==REAL(0.0) ? REAL(1.0) : REAL(3.0)*(st-t*ct)/(t*t*t) );
    const real s = (sld - solvent_sld) * form_volume(req_minor, req_major, rpolar);

    return REAL(1.0e-4) * fq * fq * s * s;
}

