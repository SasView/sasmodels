double
sphere_form(double q, double radius, double sld, double solvent_sld)
{
    const double qr = q*radius;
    double sn, cn;
    SINCOS(qr, sn, cn);
    const double vol = 1.333333333333333*M_PI*radius*radius*radius;
    const double bes = (qr == 0.0 ? 1.0 : 3.0*(sn-qr*cn)/(qr*qr*qr));
    const double fq = bes * (sld - solvent_sld)*vol;
    return 1.0e-4*fq*fq;
}

// Do we even need the paracrystal form?
// The only difference is the return value at r=0
double
sphere_form_paracrystal(double q, double radius, double delrho)
{
    const double qr = q*radius;
    const double vol = 1.333333333333333*M_PI*radius*radius*radius;
	const double f = vol*delrho*delrho*1.0e-4;
    double sn, cn;
    SINCOS(qr, sn, cn);
	const double bes = (qr == 0.0 ? f : 3.0*(sn-qr*cn)/(qr*qr*qr));
	const double fq = bes*delrho*vol;
	return fq*fq * 1.0e-4;
}