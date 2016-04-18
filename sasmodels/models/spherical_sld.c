//Headers
static double form_volume(double thick_inter[],
    double thick_flat_[],
    double core_radius,
    int n_shells);

static double Iq(double q,
    int n_shells,
    double thick_inter[],
    double func_inter[],
    double sld_core,
    double sld_solvent,
    double sld_flat[],
    double thick_flat[],
    double nu_inter[],
    int npts_inter,
    double core_radius);

static double Iqxy(double qx, double qy,
    int n_shells,
    double thick_inter[],
    double func_inter[],
    double sld_core,
    double sld_solvent,
    double sld_flat[],
    double thick_flat[],
    double nu_inter[],
    int npts_inter,
    double core_radius);

//Main code
static double form_volume(double thick_inter[],
    double thick_flat_[],
    double core_radius,
    int n)
{
    double radius = 0.0;
    int i;
    double r = core_radius;
    for (i=0; i < n; i++) {
        r += thick_inter[i];
        r += thick_flat[i];
    }
    return M_4PI_3*cube(r);
}


static double sphere_sld_kernel(double dp[], double q) {
  int n = dp[0];
  int i,j,k;

  double scale = dp[1];
  double thick_inter_core = dp[2];
  double sld_core = dp[4];
  double sld_solv = dp[5];
  double background = dp[6];
  double npts = dp[57]; //number of sub_layers in each interface
  double nsl=npts;//21.0; //nsl = Num_sub_layer:  must be ODD double number
  int n_s;

  double sld_i,sld_f,dz,bes,fun,f,vol,qr,r,contr,f2;
  double sign,slope=0.0;
  double pi;

  double total_thick=0.0;

  int fun_type[12];
  double sld[12];
  double thick_inter[12];
  double thick[12];
  double fun_coef[12];

  fun_type[0] = dp[3];
  fun_coef[0] = fabs(dp[58]);
  for (i =1; i<=n; i++){
    sld[i] = dp[i+6];
    thick_inter[i]= dp[i+16];
    thick[i] = dp[i+26];
    fun_type[i] = dp[i+36];
    fun_coef[i] = fabs(dp[i+46]);
    total_thick += thick[i];
    total_thick += thick_inter[i];
  }
  sld[0] = sld_core;
  sld[n+1] = sld_solv;
  thick[0] = dp[59];
  thick[n+1] = total_thick/5.0;
  thick_inter[0] = thick_inter_core;
  thick_inter[n+1] = 0.0;
  fun_coef[n+1] = 0.0;
  pi = 4.0*atan(1.0);
  f = 0.0;
  r = 0.0;
  vol = 0.0;
  //vol_pre = 0.0;
  //vol_sub = 0.0;
  sld_f = sld_core;

  //floor_nsl = floor(nsl/2.0);

  dz = 0.0;
  // iteration for # of shells + core + solvent
  for (i=0;i<=n+1; i++){
    //iteration for N sub-layers
    //if (fabs(thick[i]) <= 1e-24){
    //   continue;
    //}
    // iteration for flat and interface
    for (j=0;j<2;j++){
      // iteration for sub_shells in the interface
      // starts from #1 sub-layer
      for (n_s=1;n_s<=nsl; n_s++){
        // for solvent, it doesn't have an interface
        if (i==n+1 && j==1)
          break;
        // for flat layers
        if (j==0){
          dz = thick[i];
          sld_i = sld[i];
          slope = 0.0;
        }
        // for interfacial sub_shells
        else{
          dz = thick_inter[i]/nsl;
          // find sld_i at the outer boundary of sub-layer #n_s
          sld_i = intersldfunc(fun_type[i],nsl, n_s, fun_coef[i], sld[i], sld[i+1]);
          // calculate slope
          slope= (sld_i -sld_f)/dz;
        }
        contr = sld_f-slope*r;
        // iteration for the left and right boundary of the shells(or sub_shells)
        for (k=0; k<2; k++){
          // At r=0, the contribution to I is zero, so skip it.
          if ( i == 0 && j == 0 && k == 0){
            continue;
          }
          // On the top of slovent there is no interface; skip it.
          if (i == n+1 && k == 1){
            continue;
          }
          // At the right side (outer) boundary
          if ( k == 1){
            sign = 1.0;
            r += dz;
          }
          // At the left side (inner) boundary
          else{
            sign = -1.0;
          }
          qr = q * r;
          fun = 0.0;

          if(qr == 0.0){
            // sigular point
            bes = sign * 1.0;
          }
          else{
            // for flat sub-layer
            //TODO: Single precision calculation most likely fails here
            //bes = sign *  3.0 * (sin(qr) - qr * cos(qr)) / (qr * qr * qr);
            bes = sign *  sph_j1c(qr);
            if (fabs(slope) > 0.0 ){
              //fun = sign * 3.0 * r * (2.0*qr*sin(qr)-((qr*qr)-2.0)*cos(qr))/(qr * qr * qr * qr);
              fun = sign * r * sph_j1c(qr) + sign * 3.0 * sin(qr)/(qr * qr * q )
                + sign * 6.0 * cos(qr)/(qr * qr * qr * q);
            }
          }

          //Some initial optimization tries
          /*bes = (qr == 0.0 ? sign * 1.0 : sign *  3.0 * (sin(qr) - qr * cos(qr)) / (qr * qr * qr));
          //TODO: Will have to chnage this function
          if (qr!= 0.0 && fabs(slope) > 0.0 ){
            fun = sign * 3.0 * r * (2.0*qr*sin(qr)-((qr*qr)-2.0)*cos(qr))/(qr * qr * qr * qr);
          }*/

          // update total volume
          vol = 4.0 * pi / 3.0 * r * r * r;
          // we won't do the following volume correction for now.
          // substrate empty area of volume
          //if (k == 1 && fabs(sld_in[i]-sld_solv) < 1e-04*fabs(sld_solv) && fun_type[i]==0){
          //  vol_sub += (vol_pre - vol);
          //}
          f += vol * (bes * contr + fun * slope);
        }
        // remember this sld as sld_f
        sld_f = sld_i;
        // no sub-layer iteration (n_s loop) for the flat layer
        if (j==0)
          break;
      }
    }
  }
  //vol += vol_sub;
  f2 = f * f / vol;
  //f2 *= scale;
  //f2 += background;
  //free(fun_type);
  //free(sld);
  //free(thick_inter);
  //free(thick);
  //free(fun_coef);

  return (f2);
}


/**
 * Function to evaluate 1D SphereSLD function
 * @param q: q-value
 * @return: function value
 */
double Iq(double q,
    int n_shells,
    double thick_inter[],
    double func_inter[],
    double sld_core,
    double sld_solvent,
    double sld_flat[],
    double thick_flat[],
    double nu_inter[],
    int npts_inter,
    double core_radius
    ) {

    //printf("Number of points %d\n",npts_inter);
    double intensity;
    //TODO: Remove this container at later stage. It is only kept to minimize stupid errors now
    double dp[60];
    dp[0] = n_shells;
    //This is scale will also have to be removed at some stage
    dp[1] = 1.0;
    dp[2] = thick_inter_0;
    dp[3] = func_inter_0;
    dp[4] = sld_core;
    dp[5] = sld_solvent;
    dp[6] = 0.0;

    for (i=0; i<n; i++){
        dp[i+7] = sld_flat[i];
        dp[i+17] = thick_inter[i];
        dp[i+27] = thick_flat[i];
        dp[i+37] = func_inter[i];
        dp[i+47] = nu_inter[i];
    }

    dp[57] = npts_inter;
    dp[58] = nu_inter_0;
    dp[59] = rad_core_0;

    intensity = 1.0e-4*sphere_sld_kernel(dp,q);
    //printf("%10d\n",intensity);
    return intensity;
}

/**
 * Function to evaluate 2D SphereSLD function
 * @param q_x: value of Q along x
 * @param q_y: value of Q along y
 * @return: function value
 */
double Iqxy(double qx, double qy,
    int n_shells,
    double thick_inter[],
    double func_inter[],
    double sld_core,
    double sld_solvent,
    double sld_flat[],
    double thick_flat[],
    double nu_inter[],
    int npts_inter,
    double core_radius
    ) {

    double q = sqrt(qx*qx + qy*qy);
    return Iq(q, n_shells, thick_inter[], func_inter[], sld_core, sld_solvent,
    sld_flat[], thick_flat[],nu_inter[], npts_inter, core_radius)

}

