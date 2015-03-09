/*	LamellarCaille kernel - allows for name changes of passed parameters ...

*/

double Iq(double qval,
      double del,
      double Nlayers, 
      double dd,
	  double Cp, 
      double sld,
      double solvent_sld);

double Iq(double qval,
      double del,
      double Nlayers, 
      double dd,
	  double Cp, 
      double sld,
      double solvent_sld)
{
  double contr,NN;   //local variables of coefficient wave
  double inten,Pq,Sq,alpha,temp,t2;
  //double dQ, dQDefault, t1, t3;
  int ii,NNint;
  // from wikipedia 0.577215664901532860606512090082402431042159335
  const double Euler = 0.577215664901533;   // Euler's constant, increased sig figs for new models Feb 2015
  //dQDefault = 0.0;    //[=] 1/A, q-resolution, default value
  //dQ = dQDefault; // REMOVED UNUSED dQ calculations for new models Feb 2015

  NN = trunc(Nlayers);    //be sure that NN is an integer
  
  contr = sld - solvent_sld;

  Pq = 2.0*contr*contr/qval/qval*(1.0-cos(qval*del));

  NNint = (int)NN;    //cast to an integer for the loop
  ii=0;
  Sq = 0.0;
  // the vital "=" in ii<=  added March 2015
  for(ii=1;ii<=(NNint-1);ii+=1) {

    //fii = (double)ii;   //do I really need to do this? - unused variable, removed 18Feb2015

    temp = 0.0;
    alpha = Cp/4.0/M_PI/M_PI*(log(M_PI*ii) + Euler);
    //t1 = 2.0*dQ*dQ*dd*dd*alpha;
    t2 = 2.0*qval*qval*dd*dd*alpha;
    //t3 = dQ*dQ*dd*dd*ii*ii;

    temp = 1.0-ii/NN;
    //temp *= cos(dd*qval*ii/(1.0+t1));
    temp *= cos(dd*qval*ii);
    //temp *= exp(-1.0*(t2 + t3)/(2.0*(1.0+t1)) );
    temp *= exp(-t2/2.0 );
    //temp /= sqrt(1.0+t1);

    Sq += temp;
  }

  Sq *= 2.0;
  Sq += 1.0;

  inten = 2.0*M_PI*Pq*Sq/(dd*qval*qval);

  inten *= 1.0e-04;   // 1/A to 1/cm

  return(inten);
}

