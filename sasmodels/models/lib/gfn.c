//
//     FUNCTION gfn4:    CONTAINS F(Q,A,B,MU)**2  AS GIVEN
//                       BY (53) & (58-59) IN CHEN AND
//                       KOTLARCHYK REFERENCE
//
//       <OBLATE ELLIPSOID>
// function gfn4 for oblate ellipsoids
static double
gfn4(double xx, double crmaj, double crmin, double trmaj, double trmin, double delpc, double delps, double qq)
{
	// local variables
	double aa,bb,u2,ut2,uq,ut,vc,vt,siq,sit,gfnc,gfnt,tgfn,gfn4,pi43,Pi;

	Pi = 4.0*atan(1.0);
	pi43=4.0/3.0*Pi;
  	aa = crmaj;
 	bb = crmin;
 	u2 = (bb*bb*xx*xx + aa*aa*(1.0-xx*xx));
 	ut2 = (trmin*trmin*xx*xx + trmaj*trmaj*(1.0-xx*xx));
   	uq = sqrt(u2)*qq;
 	ut= sqrt(ut2)*qq;
	vc = pi43*aa*aa*bb;
   	vt = pi43*trmaj*trmaj*trmin;
   	if (uq == 0.0){
   		siq = 1.0/3.0;
   	}else{
   		siq = (sin(uq)/uq/uq - cos(uq)/uq)/uq;
   	}
   	if (ut == 0.0){
   		sit = 1.0/3.0;
   	}else{
   		sit = (sin(ut)/ut/ut - cos(ut)/ut)/ut;
   	}
   	gfnc = 3.0*siq*vc*delpc;
  	gfnt = 3.0*sit*vt*delps;
  	tgfn = gfnc+gfnt;
  	gfn4 = tgfn*tgfn;

  	return (gfn4);
}
