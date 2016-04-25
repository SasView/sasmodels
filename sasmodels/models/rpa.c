double Iq(double q, double case_num,
    double Na, double Phia, double va, double a_sld, double ba,
    double Nb, double Phib, double vb, double b_sld, double bb,
    double Nc, double Phic, double vc, double c_sld, double bc,
    double Nd, double Phid, double vd, double d_sld, double bd,
    double Kab, double Kac, double Kad,
    double Kbc, double Kbd, double Kcd
    );

double Iqxy(double qx, double qy, double case_num,
    double Na, double Phia, double va, double a_sld, double ba,
    double Nb, double Phib, double vb, double b_sld, double bb,
    double Nc, double Phic, double vc, double c_sld, double bc,
    double Nd, double Phid, double vd, double d_sld, double bd,
    double Kab, double Kac, double Kad,
    double Kbc, double Kbd, double Kcd
    );

double form_volume(void);

double form_volume(void)
{
    return 1.0;
}

double Iq(double q, double case_num,
    double Na, double Phia, double va, double La, double ba,
    double Nb, double Phib, double vb, double Lb, double bb,
    double Nc, double Phic, double vc, double Lc, double bc,
    double Nd, double Phid, double vd, double Ld, double bd,
    double Kab, double Kac, double Kad,
    double Kbc, double Kbd, double Kcd
    ) {
  int icase = (int)case_num;

#if 0  // Sasview defaults
  if (icase <= 1) {
    Na=Nb=1000.0;
    Phia=Phib=0.0000001;
    Kab=Kac=Kad=Kbc=Kbd=-0.0004;
    La=Lb=1.0e-12;
    va=vb=100.0;
    ba=bb=5.0;
  } else if (icase <= 4) {
    Phia=0.0000001;
    Kab=Kac=Kad=-0.0004;
    La=1.0e-12;
    va=100.0;
    ba=5.0;
  }
#else
  if (icase <= 1) {
    Na=Nb=0.0;
    Phia=Phib=0.0;
    Kab=Kac=Kad=Kbc=Kbd=0.0;
    La=Lb=Ld;
    va=vb=vd;
    ba=bb=0.0;
  } else if (icase <= 4) {
    Na = 0.0;
    Phia=0.0;
    Kab=Kac=Kad=0.0;
    La=Ld;
    va=vd;
    ba=0.0;
  }
#endif

  const double Xa = q*q*ba*ba*Na/6.0;
  const double Xb = q*q*bb*bb*Nb/6.0;
  const double Xc = q*q*bc*bc*Nc/6.0;
  const double Xd = q*q*bd*bd*Nd/6.0;

  // limit as Xa goes to 0 is 1
  const double Pa = Xa==0 ? 1.0 : -expm1(-Xa)/Xa;
  const double Pb = Xb==0 ? 1.0 : -expm1(-Xb)/Xb;
  const double Pc = Xc==0 ? 1.0 : -expm1(-Xc)/Xc;
  const double Pd = Xd==0 ? 1.0 : -expm1(-Xd)/Xd;

  // limit as Xa goes to 0 is 1
  const double Paa = Xa==0 ? 1.0 : 2.0*(1.0-Pa)/Xa;
  const double Pbb = Xb==0 ? 1.0 : 2.0*(1.0-Pb)/Xb;
  const double Pcc = Xc==0 ? 1.0 : 2.0*(1.0-Pc)/Xc;
  const double Pdd = Xd==0 ? 1.0 : 2.0*(1.0-Pd)/Xd;


  // Note: S0ij only defined for copolymers; otherwise set to zero
  // 0: C/D     binary mixture
  // 1: C-D     diblock copolymer
  // 2: B/C/D   ternery mixture
  // 3: B/C-D   binary mixture,1 homopolymer, 1 diblock copolymer
  // 4: B-C-D   triblock copolymer
  // 5: A/B/C/D quaternary mixture
  // 6: A/B/C-D ternery mixture, 2 homopolymer, 1 diblock copolymer
  // 7: A/B-C-D binary mixture, 1 homopolymer, 1 triblock copolymer
  // 8: A-B/C-D binary mixture, 2 diblock copolymer
  // 9: A-B-C-D tetra-block copolymer
#if 0
  const double S0aa = icase<5
                      ? 1.0 : Na*Phia*va*Paa;
  const double S0bb = icase<2
                      ? 1.0 : Nb*Phib*vb*Pbb;
  const double S0cc = Nc*Phic*vc*Pcc;
  const double S0dd = Nd*Phid*vd*Pdd;
  const double S0ab = icase<8
                      ? 0.0 : sqrt(Na*va*Phia*Nb*vb*Phib)*Pa*Pb;
  const double S0ac = icase<9
                      ? 0.0 : sqrt(Na*va*Phia*Nc*vc*Phic)*Pa*Pc*exp(-Xb);
  const double S0ad = icase<9
                      ? 0.0 : sqrt(Na*va*Phia*Nd*vd*Phid)*Pa*Pd*exp(-Xb-Xc);
  const double S0bc = (icase!=4 && icase!=7 && icase!= 9)
                      ? 0.0 : sqrt(Nb*vb*Phib*Nc*vc*Phic)*Pb*Pc;
  const double S0bd = (icase!=4 && icase!=7 && icase!= 9)
                      ? 0.0 : sqrt(Nb*vb*Phib*Nd*vd*Phid)*Pb*Pd*exp(-Xc);
  const double S0cd = (icase==0 || icase==2 || icase==5)
                      ? 0.0 : sqrt(Nc*vc*Phic*Nd*vd*Phid)*Pc*Pd;
#else  // sasview equivalent
//printf("Xc=%g, S0cc=%g*%g*%g*%g\n",Xc,Nc,Phic,vc,Pcc);
  double S0aa = Na*Phia*va*Paa;
  double S0bb = Nb*Phib*vb*Pbb;
  double S0cc = Nc*Phic*vc*Pcc;
  double S0dd = Nd*Phid*vd*Pdd;
  double S0ab = sqrt(Na*va*Phia*Nb*vb*Phib)*Pa*Pb;
  double S0ac = sqrt(Na*va*Phia*Nc*vc*Phic)*Pa*Pc*exp(-Xb);
  double S0ad = sqrt(Na*va*Phia*Nd*vd*Phid)*Pa*Pd*exp(-Xb-Xc);
  double S0bc = sqrt(Nb*vb*Phib*Nc*vc*Phic)*Pb*Pc;
  double S0bd = sqrt(Nb*vb*Phib*Nd*vd*Phid)*Pb*Pd*exp(-Xc);
  double S0cd = sqrt(Nc*vc*Phic*Nd*vd*Phid)*Pc*Pd;
switch(icase){
  case 0:
    S0aa=0.000001;
    S0ab=0.000002;
    S0ac=0.000003;
    S0ad=0.000004;
    S0bb=0.000005;
    S0bc=0.000006;
    S0bd=0.000007;
    S0cd=0.000008;
    break;
  case 1:
    S0aa=0.000001;
    S0ab=0.000002;
    S0ac=0.000003;
    S0ad=0.000004;
    S0bb=0.000005;
    S0bc=0.000006;
    S0bd=0.000007;
    break;
  case 2:
    S0aa=0.000001;
    S0ab=0.000002;
    S0ac=0.000003;
    S0ad=0.000004;
    S0bc=0.000005;
    S0bd=0.000006;
    S0cd=0.000007;
    break;
  case 3:
    S0aa=0.000001;
    S0ab=0.000002;
    S0ac=0.000003;
    S0ad=0.000004;
    S0bc=0.000005;
    S0bd=0.000006;
    break;
  case 4:
    S0aa=0.000001;
    S0ab=0.000002;
    S0ac=0.000003;
    S0ad=0.000004;
    break;
  case 5:
    S0ab=0.000001;
    S0ac=0.000002;
    S0ad=0.000003;
    S0bc=0.000004;
    S0bd=0.000005;
    S0cd=0.000006;
    break;
  case 6:
    S0ab=0.000001;
    S0ac=0.000002;
    S0ad=0.000003;
    S0bc=0.000004;
    S0bd=0.000005;
    break;
  case 7:
    S0ab=0.000001;
    S0ac=0.000002;
    S0ad=0.000003;
    break;
  case 8:
    S0ac=0.000001;
    S0ad=0.000002;
    S0bc=0.000003;
    S0bd=0.000004;
    break;
  default : //case 9:
    break;
  }
#endif

  // eq 12a: \kappa_{ij}^F = \chi_{ij}^F - \chi_{i0}^F - \chi_{j0}^F
  const double Kaa = 0.0;
  const double Kbb = 0.0;
  const double Kcc = 0.0;
  //const double Kdd = 0.0;
  const double Zaa = Kaa - Kad - Kad;
  const double Zab = Kab - Kad - Kbd;
  const double Zac = Kac - Kad - Kcd;
  const double Zbb = Kbb - Kbd - Kbd;
  const double Zbc = Kbc - Kbd - Kcd;
  const double Zcc = Kcc - Kcd - Kcd;
//printf("Za: %10.5g %10.5g %10.5g\n", Zaa, Zab, Zac);
//printf("Zb: %10.5g %10.5g %10.5g\n", Zab, Zbb, Zbc);
//printf("Zc: %10.5g %10.5g %10.5g\n", Zac, Zbc, Zcc);

  // T = inv(S0)
  const double DenT = (- S0ac*S0bb*S0ac + S0ab*S0bc*S0ac + S0ac*S0ab*S0bc
                       - S0aa*S0bc*S0bc - S0ab*S0ab*S0cc + S0aa*S0bb*S0cc);
  const double T11 = (-S0bc*S0bc + S0bb*S0cc)/DenT;
  const double T12 = ( S0ac*S0bc - S0ab*S0cc)/DenT;
  const double T13 = (-S0ac*S0bb + S0ab*S0bc)/DenT;
  const double T22 = (-S0ac*S0ac + S0aa*S0cc)/DenT;
  const double T23 = ( S0ac*S0ab - S0aa*S0bc)/DenT;
  const double T33 = (-S0ab*S0ab + S0aa*S0bb)/DenT;

//printf("T1: %10.5g %10.5g %10.5g\n", T11, T12, T13);
//printf("T2: %10.5g %10.5g %10.5g\n", T12, T22, T23);
//printf("T3: %10.5g %10.5g %10.5g\n", T13, T23, T33);

  // eq 18e: m = 1/(S0_{dd} - s0^T inv(S0) s0)
  const double ZZ = S0ad*(T11*S0ad + T12*S0bd + T13*S0cd)
                  + S0bd*(T12*S0ad + T22*S0bd + T23*S0cd)
                  + S0cd*(T13*S0ad + T23*S0bd + T33*S0cd);

  const double m=1.0/(S0dd-ZZ);

  // eq 18d: Y = inv(S0)s0 + e
  const double Y1 = T11*S0ad + T12*S0bd + T13*S0cd + 1.0;
  const double Y2 = T12*S0ad + T22*S0bd + T23*S0cd + 1.0;
  const double Y3 = T13*S0ad + T23*S0bd + T33*S0cd + 1.0;

  // N = mYY^T + \kappa^F
  const double N11 = m*Y1*Y1 + Zaa;
  const double N12 = m*Y1*Y2 + Zab;
  const double N13 = m*Y1*Y3 + Zac;
  const double N22 = m*Y2*Y2 + Zbb;
  const double N23 = m*Y2*Y3 + Zbc;
  const double N33 = m*Y3*Y3 + Zcc;

//printf("N1: %10.5g %10.5g %10.5g\n", N11, N12, N13);
//printf("N2: %10.5g %10.5g %10.5g\n", N12, N22, N23);
//printf("N3: %10.5g %10.5g %10.5g\n", N13, N23, N33);
//printf("S0a: %10.5g %10.5g %10.5g\n", S0aa, S0ab, S0ac);
//printf("S0b: %10.5g %10.5g %10.5g\n", S0ab, S0bb, S0bc);
//printf("S0c: %10.5g %10.5g %10.5g\n", S0ac, S0bc, S0cc);

  // M = I + S0 N
  const double Maa = N11*S0aa + N12*S0ab + N13*S0ac + 1.0;
  const double Mab = N11*S0ab + N12*S0bb + N13*S0bc;
  const double Mac = N11*S0ac + N12*S0bc + N13*S0cc;
  const double Mba = N12*S0aa + N22*S0ab + N23*S0ac;
  const double Mbb = N12*S0ab + N22*S0bb + N23*S0bc + 1.0;
  const double Mbc = N12*S0ac + N22*S0bc + N23*S0cc;
  const double Mca = N13*S0aa + N23*S0ab + N33*S0ac;
  const double Mcb = N13*S0ab + N23*S0bb + N33*S0bc;
  const double Mcc = N13*S0ac + N23*S0bc + N33*S0cc + 1.0;
//printf("M1: %10.5g %10.5g %10.5g\n", Maa, Mab, Mac);
//printf("M2: %10.5g %10.5g %10.5g\n", Mba, Mbb, Mbc);
//printf("M3: %10.5g %10.5g %10.5g\n", Mca, Mcb, Mcc);

  // Q = inv(M) = inv(I + S0 N)
  const double DenQ = (+ Maa*Mbb*Mcc - Maa*Mbc*Mcb - Mab*Mba*Mcc
                       + Mab*Mbc*Mca + Mac*Mba*Mcb - Mac*Mbb*Mca);

  const double Q11 = ( Mbb*Mcc - Mbc*Mcb)/DenQ;
  const double Q12 = (-Mab*Mcc + Mac*Mcb)/DenQ;
  const double Q13 = ( Mab*Mbc - Mac*Mbb)/DenQ;
  //const double Q21 = (-Mba*Mcc + Mbc*Mca)/DenQ;
  const double Q22 = ( Maa*Mcc - Mac*Mca)/DenQ;
  const double Q23 = (-Maa*Mbc + Mac*Mba)/DenQ;
  //const double Q31 = ( Mba*Mcb - Mbb*Mca)/DenQ;
  //const double Q32 = (-Maa*Mcb + Mab*Mca)/DenQ;
  const double Q33 = ( Maa*Mbb - Mab*Mba)/DenQ;

//printf("Q1: %10.5g %10.5g %10.5g\n", Q11, Q12, Q13);
//printf("Q2: %10.5g %10.5g %10.5g\n", Q21, Q22, Q23);
//printf("Q3: %10.5g %10.5g %10.5g\n", Q31, Q32, Q33);
  // eq 18c: inv(S) = inv(S0) + mYY^T + \kappa^F
  // eq A1 in the appendix
  // To solve for S, use:
  //      S = inv(inv(S^0) + N) inv(S^0) S^0
  //        = inv(S^0 inv(S^0) + N) S^0
  //        = inv(I + S^0 N) S^0
  //        = Q S^0
  const double S11 = Q11*S0aa + Q12*S0ab + Q13*S0ac;
  const double S12 = Q12*S0aa + Q22*S0ab + Q23*S0ac;
  const double S13 = Q13*S0aa + Q23*S0ab + Q33*S0ac;
  const double S22 = Q12*S0ab + Q22*S0bb + Q23*S0bc;
  const double S23 = Q13*S0ab + Q23*S0bb + Q33*S0bc;
  const double S33 = Q13*S0ac + Q23*S0bc + Q33*S0cc;
  // If the full S is needed...it isn't since Ldd = (rho_d - rho_d) = 0 below
  //const double S14=-S11-S12-S13;
  //const double S24=-S12-S22-S23;
  //const double S34=-S13-S23-S33;
  //const double S44=S11+S22+S33 + 2.0*(S12+S13+S23);

  // eq 12 of Akcasu, 1990: I(q) = L^T S L
  // Note: eliminate cases without A and B polymers by setting Lij to 0
  // Note: 1e-13 to convert from fm to cm for scattering length
  const double sqrt_Nav=sqrt(6.022045e+23) * 1.0e-13;
  const double Lad = icase<5 ? 0.0 : (La/va - Ld/vd)*sqrt_Nav;
  const double Lbd = icase<2 ? 0.0 : (Lb/vb - Ld/vd)*sqrt_Nav;
  const double Lcd = (Lc/vc - Ld/vd)*sqrt_Nav;

  const double result=Lad*Lad*S11 + Lbd*Lbd*S22 + Lcd*Lcd*S33
                    + 2.0*(Lad*Lbd*S12 + Lbd*Lcd*S23 + Lad*Lcd*S13);

  return result;

}

double Iqxy(double qx, double qy,
    double case_num,
    double Na, double Phia, double va, double a_sld, double ba,
    double Nb, double Phib, double vb, double b_sld, double bb,
    double Nc, double Phic, double vc, double c_sld, double bc,
    double Nd, double Phid, double vd, double d_sld, double bd,
    double Kab, double Kac, double Kad,
    double Kbc, double Kbd, double Kcd
    )
{
    double q = sqrt(qx*qx + qy*qy);
    return Iq(q,
        case_num,
        Na, Phia, va, a_sld, ba,
        Nb, Phib, vb, b_sld, bb,
        Nc, Phic, vc, c_sld, bc,
        Nd, Phid, vd, d_sld, bd,
        Kab, Kac, Kad,
        Kbc, Kbd, Kcd);
}