/*
The wrapper for gamma function from OpenCL and standard libraries
The OpenCL gamma function fails miserably on values lower than 1.0
while works fine on larger values.
We use gamma definition Gamma(t + 1) = t * Gamma(t) to compute
to function for values lower than 1.0. Namely Gamma(t) = 1/t * Gamma(t + 1)
*/

inline double sas_gamma( double x) { return tgamma(x+1)/x; }