/*	Lamellar_ParaCrystal - Pedersen's model

*/
double Iq(double qval,
      double th,
      double Nlayers, 
	  double davg, 
	  double pd,
      double sld,
      double solvent_sld);
double paraCryst_sn(double ww, double qval, double davg, long Nlayers, double an);
double paraCryst_an(double ww, double qval, double davg, long Nlayers);

double Iq(double qval,
      double th,
      double Nlayers, 
	  double davg, 
	  double pd,
      double sld,
      double solvent_sld)
{
    
	double inten,contr,xn;
	double xi,ww,Pbil,Znq,Snq,an;
	long n1,n2;
	
	contr = sld - solvent_sld;
	//get the fractional part of Nlayers, to determine the "mixing" of N's
	
	n1 = (long)trunc(Nlayers);		//rounds towards zero
	n2 = n1 + 1;
	xn = (double)n2 - Nlayers;			//fractional contribution of n1
	
	ww = exp(-qval*qval*pd*pd*davg*davg/2.0);
	
	//calculate the n1 contribution
	an = paraCryst_an(ww,qval,davg,n1);
	Snq = paraCryst_sn(ww,qval,davg,n1,an);
	
	Znq = xn*Snq;
	
	//calculate the n2 contribution
	an = paraCryst_an(ww,qval,davg,n2);
	Snq = paraCryst_sn(ww,qval,davg,n2,an);
	
	Znq += (1.0-xn)*Snq;
	
	//and the independent contribution
	Znq += (1.0-ww*ww)/(1.0+ww*ww-2.0*ww*cos(qval*davg));
	
	//the limit when Nlayers approaches infinity
//	Zq = (1-ww^2)/(1+ww^2-2*ww*cos(qval*davg))
	
	xi = th/2.0;		//use 1/2 the bilayer thickness
	Pbil = (sin(qval*xi)/(qval*xi))*(sin(qval*xi)/(qval*xi));
	
	inten = 2.0*M_PI*contr*contr*Pbil*Znq/(qval*qval);
	inten *= 1.0e-04;
	
	return(inten);
}

// functions for the lamellar paracrystal model
double
paraCryst_sn(double ww, double qval, double davg, long Nlayers, double an) {
	
	double Snq;
	
	Snq = an/( (double)Nlayers*pow((1.0+ww*ww-2.0*ww*cos(qval*davg)),2) );
	
	return(Snq);
}

double
paraCryst_an(double ww, double qval, double davg, long Nlayers) {
	
	double an;
	
	an = 4.0*ww*ww - 2.0*(ww*ww*ww+ww)*cos(qval*davg);
	an -= 4.0*pow(ww,(Nlayers+2))*cos((double)Nlayers*qval*davg);
	an += 2.0*pow(ww,(Nlayers+3))*cos((double)(Nlayers-1)*qval*davg);
	an += 2.0*pow(ww,(Nlayers+1))*cos((double)(Nlayers+1)*qval*davg);
	
	return(an);
}

