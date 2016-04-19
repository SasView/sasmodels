static double form_volume(
    int n_shells,
    double radius_core,
    double thick_inter0,
    double thick_flat[],
    double thick_inter[])
{
    int i;
    double r = radius_core;
    r += thick_inter0;
    for (i=0; i < n_shells; i++) {
        r += thick_inter[i];
        r += thick_flat[i];
    }
    return M_4PI_3*cube(r);
}


static double sphere_sld_kernel(
    double q,
    int n_shells,
    int npts_inter,
    double radius_core,
    double sld_core,
    double sld_solvent,
    double func_inter_core,
    double thick_inter_core,
    double nu_inter_core,
    double sld_flat[],
    double thick_flat[],
    double func_inter[],
    double thick_inter[],
    double nu_inter[] ) {

    int i,j,k;
    int n_s;

    double sld_i,sld_f,dz,bes,fun,f,vol,qr,r,contr,f2;
    double sign,slope=0.0;
    double pi;

    double total_thick=0.0;

    //TODO: This part can be further simplified
    int fun_type[12];
    double sld[12];
    double thick_internal[12];
    double thick[12];
    double fun_coef[12];

    fun_type[0] = func_inter_core;
    fun_coef[0] = fabs(nu_inter_core);
    sld[0] = sld_core;
    thick[0] = radius_core;
    thick_internal[0] = thick_inter_core;

    for (i =1; i<=n_shells; i++){
        sld[i] = sld_flat[i-1];
        thick_internal[i]= thick_inter[i-1];
        thick[i] = thick_flat[i-1];
        fun_type[i] = func_inter[i-1];
        fun_coef[i] = fabs(nu_inter[i-1]);
        total_thick += thick[i];
        total_thick += thick_internal[i]; //doesn't account for core layer
    }

    sld[n_shells+1] = sld_solvent;
    thick[n_shells+1] = total_thick/5.0;
    thick_internal[n_shells+1] = 0.0;
    fun_coef[n_shells+1] = 0.0;
    fun_type[n_shells+1] = 0;

    pi = 4.0*atan(1.0);
    f = 0.0;
    r = 0.0;
    vol = 0.0;
    sld_f = sld_core;

    dz = 0.0;
    // iteration for # of shells + core + solvent
    for (i=0;i<=n_shells+1; i++){
        // iteration for flat and interface
        for (j=0;j<2;j++){
            // iteration for sub_shells in the interface
            // starts from #1 sub-layer
            for (n_s=1;n_s<=npts_inter; n_s++){
                // for solvent, it doesn't have an interface
                if (i==n_shells+1 && j==1)
                    break;
                // for flat layers
                if (j==0){
                    dz = thick[i];
                    sld_i = sld[i];
                    slope = 0.0;
                }
                // for interfacial sub_shells
                else{
                    dz = thick_internal[i]/npts_inter;
                    // find sld_i at the outer boundary of sub-layer #n_s
                    sld_i = intersldfunc(fun_type[i], npts_inter, n_s,
                            fun_coef[i], sld[i], sld[i+1]);
                    // calculate slope
                    slope= (sld_i -sld_f)/dz;
                }
                contr = sld_f-slope*r;
                // iteration for the left and right boundary of the shells
                for (k=0; k<2; k++){
                    // At r=0, the contribution to I is zero, so skip it.
                    if ( i == 0 && j == 0 && k == 0){
                        continue;
                    }
                    // On the top of slovent there is no interface; skip it.
                    if (i == n_shells+1 && k == 1){
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
                        bes = sign * 1.0;
                    }
                    else{
                        // for flat sub-layer
                        //TODO: Single precision calculation fails here
                        bes = sign *  sph_j1c(qr);
                        if (fabs(slope) > 0.0 ){
                            const double qrsq = qr*qr;
                            double sinqr, cosqr;
                            SINCOS(qr, sinqr, cosqr);
                            fun = sign * 3.0 * r *
                            (2.0*qr*sinqr - (qrsq-2.0)*cosqr)/(qrsq * qrsq);
                            // In the onioon model Jae-He's formula is rederived
                            // and gives following:
                            //fun = 3.0 * sign * r *
                            //(2.0*cosqr + qr*sinqr)/(qrsq*qrsq);
                            //But this seems not to be working in this case...
                        }
                    }

                    // update total volume
                    vol = M_4PI_3 * cube(r);
                    f += vol * (bes * contr + fun * slope);
                }
                sld_f = sld_i;
                // no sub-layer iteration (n_s loop) for the flat layer
                if (j==0)
                    break;
            }
        }
    }

    f2 = f * f * 1.0e-4;
    return (f2);
}


/**
 * Function to evaluate 1D SphereSLD function
 * @param q: q-value
 * @return: function value
 */
static double Iq(double q,
    int n_shells,
    int npts_inter,
    double radius_core,
    double sld_core,
    double sld_solvent,
    double func_inter0,
    double thick_inter0,
    double nu_inter0,
    double sld_flat[],
    double thick_flat[],
    double func_inter[],
    double thick_inter[],
    double nu_inter[] ) {

    //printf("Number of points %d\n",npts_inter);
    double intensity;
    //TODO: Remove this container at later stage.
    /*double dp[60];
    dp[0] = n_shells;
    //This is scale will also have to be removed at some stage
    dp[1] = 1.0;
    dp[2] = thick_inter0;
    dp[3] = func_inter0;
    dp[4] = sld_core;
    dp[5] = sld_solvent;
    dp[6] = 0.0;

    for (int i=0; i<n_shells; i++){
        dp[i+7] = sld_flat[i];
        dp[i+17] = thick_inter[i];
        dp[i+27] = thick_flat[i];
        dp[i+37] = func_inter[i];
        dp[i+47] = nu_inter[i];
    }

    dp[57] = npts_inter;
    dp[58] = nu_inter0;
    dp[59] = radius_core;
    */
    intensity = sphere_sld_kernel(q, n_shells, npts_inter, radius_core,
                sld_core, sld_solvent, func_inter0, thick_inter0, nu_inter0,
                sld_flat, thick_flat, func_inter, thick_inter, nu_inter);
    //intensity *=1.0e-4;
    //printf("%10d\n",intensity);
    return intensity;
}

/**
 * Function to evaluate 2D SphereSLD function
 * @param q_x: value of Q along x
 * @param q_y: value of Q along y
 * @return: function value
 */

/*static double Iqxy(double qx, double qy,
    int n_shells,
    int npts_inter,
    double radius_core
    double sld_core,
    double sld_solvent,
    double sld_flat[],
    double thick_flat[],
    double func_inter[],
    double thick_inter[],
    double nu_inter[],
    ) {

    double q = sqrt(qx*qx + qy*qy);
    return Iq(q, n_shells, npts_inter, radius_core, sld_core, sld_solvent,
    sld_flat[], thick_flat[], func_inter[], thick_inter[], nu_inter[])
}*/

