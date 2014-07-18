__kernel void EllipsoidKernel(__global const real *qx, __global const real *qy, __global const int *place, __global const real *array,
__global real *final, const real scale, const real sub, const int length, const int size)//, __local real *poly)
{
    //__local real rada, radaw, radb, radbw, th, thw, ph, phw;

    __local real rada;
    int uno=0; int dos=1; int count = 0;
    int l = get_local_id(0);
    for(int j = 0; j < 4; j++)
    {
        for(int i = place[uno]; i < place[dos]; i++)
        {
            rada[count] = array[i]
            radaw[count] = array[]
        }
    uno+=2;
    dos+=2;
    count = 0;
    }


    int i = get_global_id(0);
    if(i < length){
    //do local mem?
    real sum=0.0; real vol=0.0; real norm_vol=0.0; real norm=0.0;

    for(int a = 0; a < place[0]; a++) //radius_a values
    {
        real radius_a = array[a];
        real radius_a_weight = array[a+place[0]];
        for(int b = place[1]; b < place[2]; b++) //radius_b values
        {
            real radius_b = array[b];
            real radius_b_weight = array[b+place[2]];
            for(int t = place[3]; t < place[4]; t++) //Axis_theta values
            {
                real axis_theta = array[t];
                real axis_theta_weight = array[t+place[4]];
                for(int p = place[5]; p < place[6]; p++) //axis_phi values
                {
                    real axis_phi = array[p];
                    real axis_phi_weight = array[p+place[6]];

                    real ret = 0;
                    real q = sqrt(qx[i]*qx[i] + qy[i]*qy[i]);
                    real pi = 4.0*atan(1.0);
                    real theta = axis_theta*pi/180.0;
                    real cos_val = cos(theta)*cos(axis_phi*pi/180.0)*(qx[i]/q) + sin(theta)*(qy[i]/q);

                    real arg = q*radius_b*sqrt(1.0+(cos_val*cos_val*(((radius_a*radius_a/(radius_b*radius_b))-1.0))));
                    if(arg == 0.0){
                     ret = 1.0/3.0;
                    }
                    else{
                        ret = (sin(arg)-arg*cos(arg))/(arg*arg*arg);
                    }
                    ret*=ret*9.0*sub*sub*4.0/3.0*acos(-1.0)*radius_b*radius_b*radius_a*scale*(1.0e8);
                    real _ptvalue = radius_a_weight*radius_b_weight*axis_theta_weight*radius_a*axis_phi_weight*ret*pow(radius_b, 2);
                    if(size > 1){
                        _ptvalue *= fabs(cos(axis_theta*pi/180.0));
                    }

                    sum += _ptvalue;
                    vol += radius_a_weight*radius_b_weight*pow(radius_b, 2)*radius_a;
                    norm_vol += radius_a_weight*radius_b_weight;
                    norm += radius_a_weight*radius_b_weight*axis_theta_weight*axis_phi_weight;
                }
            }
        }
    }

    if(size > 1){
        norm /= asin(1.0);
    }
    if(vol != 0.0 && norm_vol != 0.0){
        sum *= norm_vol/vol;
    }

    final[i] = sum/norm;
}
}