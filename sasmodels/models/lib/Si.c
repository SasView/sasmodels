double Si(double x);
// integral of sin(x)/x: approximated to w/i 1%
double Si(double x)
{
	int i;
	int nmax=6;
	double out;
	long power;
	double pi = 4.0*atan(1.0);

	if (x >= pi*6.2/4.0){
		double out_sin = 0.0;
		double out_cos = 0.0;
		out = pi/2.0;

		for (i=0; i<nmax-2; i+=1) {
			out_cos += pow(-1.0, i) * (double)factorial(2*i) / pow(x, 2*i+1);
			out_sin += pow(-1.0, i) * (double)factorial(2*i+1) / pow(x, 2*i+2);
		}

		out -= cos(x) * out_cos;
		out -= sin(x) * out_sin;
		return out;
	}

	out = 0.0;

	for (i=0; i<nmax; i+=1)	{
		if (i==0) {
			out += x;
			continue;
		}

		power = pow(x,(2 * i + 1));
		out += (double)pow(-1, i) * power / ((2.0 * (double)i + 1.0) * (double)factorial(2 * i + 1));

		//printf ("Si=%g %g %d\n", x, out, i);
	}

	return out;
}

int factorial(int f);
int factorial(int f)
{
    if ( f == 0 ) 
        return 1;
    return(f * factorial(f - 1));
}