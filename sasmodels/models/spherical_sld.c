//Headers
static double form_volume(double thick_inter[],
    double thick_flat_[],
    double core_radius,
    int n_shells);

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
    double core_radius);

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
  double nsl=npts;//21.0; //nsl = Num_sub_layer:  MUST ODD number in double //no other number works now
  int n_s;

  double sld_i,sld_f,dz,bes,fun,f,vol,qr,r,contr,f2;
  double sign,slope=0.0;
  double pi;

  //int* fun_type;
  //double* sld;
  //double* thick_inter;
  //double* thick;
  //double* fun_coef;

  double total_thick=0.0;

  //fun_type = (int*)malloc((n+2)*sizeof(int));
  //sld = (double*)malloc((n+2)*sizeof(double));
  //thick_inter = (double*)malloc((n+2)*sizeof(double));
  //thick = (double*)malloc((n+2)*sizeof(double));
  //fun_coef = (double*)malloc((n+2)*sizeof(double));

  //TODO: Solution to avoid mallocs but probablyu can be done better
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
              fun = sign * r * sph_j1c(qr)  +  sign * 3.0 * sin(qr)/(qr * qr * q ) + sign * 6.0 * cos(qr)/(qr * qr * qr * q);
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
    int n_shells, double thick_inter_0, int func_inter_0, double core0_sld, double solvent_sld,
    double flat1_sld, double flat2_sld, double flat3_sld, double flat4_sld, double flat5_sld,
    double flat6_sld, double flat7_sld, double flat8_sld, double flat9_sld, double flat10_sld,
    double thick_inter_1, double thick_inter_2, double thick_inter_3, double thick_inter_4, double thick_inter_5,
    double thick_inter_6, double thick_inter_7, double thick_inter_8, double thick_inter_9, double thick_inter_10,
    double thick_flat_1, double thick_flat_2, double thick_flat_3, double thick_flat_4, double thick_flat_5,
    double thick_flat_6, double thick_flat_7, double thick_flat_8, double thick_flat_9, double thick_flat_10,
    int func_inter_1, int func_inter_2, int func_inter_3, int func_inter_4, int func_inter_5,
    int func_inter_6, int func_inter_7, int func_inter_8, int func_inter_9, int func_inter_10,
    double nu_inter_1, double nu_inter_2,double nu_inter_3, double nu_inter_4, double nu_inter_5,
    double nu_inter_6, double nu_inter_7, double nu_inter_8, double nu_inter_9, double nu_inter_10,
    int npts_inter, double nu_inter_0, double rad_core_0) {

    //printf("Number of points %d\n",npts_inter);
    double intensity;
    //TODO: Remove this container at later stage. It is only kept to minimize stupid errors now
    double dp[60];
    dp[0] = n_shells;
    //This is scale will also have to be removed at some stage
    dp[1] = 1.0;
    dp[2] = thick_inter_0;
    dp[3] = func_inter_0;
    dp[4] = core0_sld;
    dp[5] = solvent_sld;
    dp[6] = 0.0;

    dp[7] = flat1_sld;
    dp[8] = flat2_sld;
    dp[9] = flat3_sld;
    dp[10] = flat4_sld;
    dp[11] = flat5_sld;
    dp[12] = flat6_sld;
    dp[13] = flat7_sld;
    dp[14] = flat8_sld;
    dp[15] = flat9_sld;
    dp[16] = flat10_sld;

    dp[17] = thick_inter_1;
    dp[18] = thick_inter_2;
    dp[19] = thick_inter_3;
    dp[20] = thick_inter_4;
    dp[21] = thick_inter_5;
    dp[22] = thick_inter_6;
    dp[23] = thick_inter_7;
    dp[24] = thick_inter_8;
    dp[25] = thick_inter_9;
    dp[26] = thick_inter_10;

    dp[27] = thick_flat_1;
    dp[28] = thick_flat_2;
    dp[29] = thick_flat_3;
    dp[30] = thick_flat_4;
    dp[31] = thick_flat_5;
    dp[32] = thick_flat_6;
    dp[33] = thick_flat_7;
    dp[34] = thick_flat_8;
    dp[35] = thick_flat_9;
    dp[36] = thick_flat_10;

    dp[37] = func_inter_1;
    dp[38] = func_inter_2;
    dp[39] = func_inter_3;
    dp[40] = func_inter_4;
    dp[41] = func_inter_5;
    dp[42] = func_inter_6;
    dp[43] = func_inter_7;
    dp[44] = func_inter_8;
    dp[45] = func_inter_9;
    dp[46] = func_inter_10;

    dp[47] = nu_inter_1;
    dp[48] = nu_inter_2;
    dp[49] = nu_inter_3;
    dp[50] = nu_inter_4;
    dp[51] = nu_inter_5;
    dp[52] = nu_inter_6;
    dp[53] = nu_inter_7;
    dp[54] = nu_inter_8;
    dp[55] = nu_inter_9;
    dp[56] = nu_inter_10;

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
    int n_shells, double thick_inter_0, int func_inter_0, double core0_sld, double solvent_sld,
    double flat1_sld, double flat2_sld, double flat3_sld, double flat4_sld, double flat5_sld,
    double flat6_sld, double flat7_sld, double flat8_sld, double flat9_sld, double flat10_sld,
    double thick_inter_1, double thick_inter_2, double thick_inter_3, double thick_inter_4, double thick_inter_5,
    double thick_inter_6, double thick_inter_7, double thick_inter_8, double thick_inter_9, double thick_inter_10,
    double thick_flat_1, double thick_flat_2, double thick_flat_3, double thick_flat_4, double thick_flat_5,
    double thick_flat_6, double thick_flat_7, double thick_flat_8, double thick_flat_9, double thick_flat_10,
    int func_inter_1, int func_inter_2, int func_inter_3, int func_inter_4, int func_inter_5,
    int func_inter_6, int func_inter_7, int func_inter_8, int func_inter_9, int func_inter_10,
    double nu_inter_1, double nu_inter_2,double nu_inter_3, double nu_inter_4, double nu_inter_5,
    double nu_inter_6, double nu_inter_7, double nu_inter_8, double nu_inter_9, double nu_inter_10,
    int npts_inter, double nu_inter_0, double rad_core_0) {

    double q = sqrt(qx*qx + qy*qy);
    return Iq(q, n_shells, thick_inter_0, func_inter_0, core0_sld, solvent_sld,
    flat1_sld, flat2_sld, flat3_sld, flat4_sld, flat5_sld,
    flat6_sld, flat7_sld, flat8_sld, flat9_sld, flat10_sld,
    thick_inter_1, thick_inter_2, thick_inter_3, thick_inter_4, thick_inter_5,
    thick_inter_6, thick_inter_7, thick_inter_8, thick_inter_9, thick_inter_10,
    thick_flat_1, thick_flat_2, thick_flat_3, thick_flat_4, thick_flat_5,
    thick_flat_6, thick_flat_7, thick_flat_8, thick_flat_9, thick_flat_10,
    func_inter_1, func_inter_2, func_inter_3, func_inter_4, func_inter_5,
    func_inter_6, func_inter_7, func_inter_8, func_inter_9, func_inter_10,
    nu_inter_1, nu_inter_2, nu_inter_3, nu_inter_4, nu_inter_5,
    nu_inter_6, nu_inter_7, nu_inter_8, nu_inter_9, nu_inter_10,
    npts_inter, nu_inter_0, rad_core_0);

    //TODO: Check if evalute rphi is not needed?

}

