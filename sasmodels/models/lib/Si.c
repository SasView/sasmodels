double Si(double x);

// integral of sin(x)/x Taylor series approximated to w/i 0.1%
double Si(double x)
{
	double out;
	double pi = 4.0*atan(1.0);

	if (x >= pi*6.2/4.0){
		double out_sin = 0.0;
		double out_cos = 0.0;
		out = pi/2.0;

		// Explicitly writing factorial values triples the speed of the calculation
		out_cos = 1/x - 2/pow(x,3) + 24/pow(x,5) - 720/pow(x,7) + 40320/pow(x,9);
		out_sin = 1/x - 6/pow(x,4) + 120/pow(x,6) - 5040/pow(x,8) + 362880/pow(x,10);

		out -= cos(x) * out_cos;
		out -= sin(x) * out_sin;
		return out;
	}

	// Explicitly writing factorial values triples the speed of the calculation
	out = x - pow(x, 3)/18 + pow(x,5)/600 - pow(x,7)/35280 + pow(x,9)/3265920;

	//printf ("Si=%g %g\n", x, out);
	return out;
}