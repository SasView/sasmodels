real ConvLens_kernel(real len, real rad, real endRad, real x, real tt, real theta)
{
    real be = 0;
	real hDist = -1.0*sqrt(fabs(endRad*endRad-rad*rad));
	real arg1 = x*cos(theta)*(endRad*tt+hDist+len/2.0);
	real arg2 = x*endRad*sin(theta)*sqrt(1.0-tt*tt);

	if(arg2 == 0) {be = 0.5;}
	else {
		be = NR_BessJ1(arg2)/arg2;
	}
	real val = cos(arg1)*(1.0-tt*tt)*be;

	return(val);
}