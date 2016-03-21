FourShell(double dp[], double q)
{
    // variables are:
    //[0] scale factor
    //[1] radius of core [�]
    //[2] SLD of the core   [�-2]
    //[3] thickness of shell 1 [�]
    //[4] SLD of shell 1
    //[5] thickness of shell 2 [�]
    //[6] SLD of shell 2
    //[7] thickness of shell 3
    //[8] SLD of shell 3
    //[9] thickness of shell 3
    //[10] SLD of shell 3
    //[11] SLD of solvent
    //[12] background   [cm-1]
    
    double x,pi;
    double scale,rcore,thick1,rhocore,rhoshel1,rhosolv,bkg;     //my local names
    double bes,f,vol,qr,contr,f2;
    double rhoshel2,thick2,rhoshel3,thick3,rhoshel4,thick4;
    
    pi = 4.0*atan(1.0);
    x=q;
    
    scale = dp[0];
    rcore = dp[1];
    rhocore = dp[2];
    thick1 = dp[3];
    rhoshel1 = dp[4];
    thick2 = dp[5];
    rhoshel2 = dp[6];   
    thick3 = dp[7];
    rhoshel3 = dp[8];
    thick4 = dp[9];
    rhoshel4 = dp[10];  
    rhosolv = dp[11];
    bkg = dp[12];
    
        // core first, then add in shells
    qr=x*rcore;
    contr = rhocore-rhoshel1;
    if(qr == 0){
        bes = 1.0;
    }else{
        bes = 3.0*(sin(qr)-qr*cos(qr))/(qr*qr*qr);
    }
    vol = 4.0*pi/3.0*rcore*rcore*rcore;
    f = vol*bes*contr;
    //now the shell (1)
    qr=x*(rcore+thick1);
    contr = rhoshel1-rhoshel2;
    if(qr == 0){
        bes = 1.0;
    }else{
        bes = 3.0*(sin(qr)-qr*cos(qr))/(qr*qr*qr);
    }
    vol = 4.0*pi/3.0*(rcore+thick1)*(rcore+thick1)*(rcore+thick1);
    f += vol*bes*contr;
    //now the shell (2)
    qr=x*(rcore+thick1+thick2);
    contr = rhoshel2-rhoshel3;
    if(qr == 0){
        bes = 1.0;
    }else{
        bes = 3.0*(sin(qr)-qr*cos(qr))/(qr*qr*qr);
    }
    vol = 4.0*pi/3.0*(rcore+thick1+thick2)*(rcore+thick1+thick2)*(rcore+thick1+thick2);
    f += vol*bes*contr;
    //now the shell (3)
    qr=x*(rcore+thick1+thick2+thick3);
    contr = rhoshel3-rhoshel4;
    if(qr == 0){
        bes = 1.0;
    }else{
        bes = 3.0*(sin(qr)-qr*cos(qr))/(qr*qr*qr);
    }
    vol = 4.0*pi/3.0*(rcore+thick1+thick2+thick3)*(rcore+thick1+thick2+thick3)*(rcore+thick1+thick2+thick3);
    f += vol*bes*contr;
    //now the shell (4)
    qr=x*(rcore+thick1+thick2+thick3+thick4);
    contr = rhoshel4-rhosolv;
    if(qr == 0){
        bes = 1.0;
    }else{
        bes = 3.0*(sin(qr)-qr*cos(qr))/(qr*qr*qr);
    }
    vol = 4.0*pi/3.0*(rcore+thick1+thick2+thick3+thick4)*(rcore+thick1+thick2+thick3+thick4)*(rcore+thick1+thick2+thick3+thick4);
    f += vol*bes*contr;
    
        
    // normalize to particle volume and rescale from [�-1] to [cm-1]
    f2 = f*f/vol*1.0e8;
    
    //scale if desired
    f2 *= scale;
    // then add in the background
    f2 += bkg;
    
    return(f2);
}

// above is cut and past with no changes from libigor four shell
static double
f_exp(double q, double r, double sld_in, double sld_out,
    double thickness, double A)
{
  const double vol = 4.0/3.0 * M_PI * r * r * r;
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
    fun = -3.0(
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
  const double vol = 4.0/3.0 * M_PI * r * r * r;
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
  const double vol = 4.0/3.0 * M_PI * r * r * r;
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
  return 4.0/3.0 * M_PI * r * r * r;
}

static double
Iq(double q, double core_sld, double core_radius, double solvent_sld,
    double n, double in_sld[], double out_sld[], double thickness[],
    double A[])
{
  int i;
  r = core_radius;
  f = f_constant(q, r, core_sld);
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
  f -= f_constant(q, r, solvent_sld);
  f2 = f * f * 1.0e-4;

  return f2;
}
