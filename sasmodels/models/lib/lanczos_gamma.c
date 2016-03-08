double lanczos_gamma(double q);
double lanczos_gamma(double q)
{
    // Lanczos approximation to the Log Gamma function.

    double x,y,tmp,ser;
    double coeff[6]=
        {76.18009172947146,     -86.50532032941677,
         24.01409824083091,     -1.231739572450155,
          0.1208650973866179e-2,-0.5395239384953e-5};

    y=x=q;
    tmp  = x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser  = 1.000000000190015;
    for (int j=0; j<=5; j++) {
        y+=1.0;
        ser += coeff[j]/y;
    }
    return -tmp+log(2.5066282746310005*ser/x);
}