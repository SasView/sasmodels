double Iq(
    double qexp, 
    double radius,
    double volumefraction,
    double k1,  
    double k2,   
    double z1,    
    double z2
    )
{
    double a;
    double b;
    double c1;
    double c2;
    double d1;
    double d2;
    int debug;
    int checkFlag;

    a = 1.0;
    b = 1.0;
    c1 = 1.0;
    d1 = 1.0;
    c2 = 1.0;
    d2 = 1.0;
    debug = 0;
    checkFlag = 0;  

    if ( z1 == z2 )
    {
 //       Y_SolveEquations(z1, k1 + k2, volumefraction, &a, &b, &c1, &d1, debug );
 //       Y_CheckSolution(z1, k1 + k2, volumefraction, &a, &b, &c1, &d1 );
 //       return SqOneYukawa(qexp, z1, k1 + k2, volumefraction, &a, &b, &c1, &d1 );
          return 0;
    }
        

// Theoretically speaking, z1 and z2 (and k1 and k2) are symetric. 
// Exchanging z1 (k1) with z2 (k2) does not altern the potential. 
// The results should be identical. However, the orginal model proposed 
// by Y. Liu treats z1 and z2 differently when implementing the numerical solution. // Hence, it is in general a good practice to require the z1 > z2. 
// The following code is added here to swap z1 (k1) with z2 (k2) if z1 < z2.

    double temp;
    if ( z1 < z2 )
    {
        temp = z1;
        z1 = z2;
        z2 = temp;
        temp = k1;
        k1 = k2;
        k2 = temp;
    }

  
  
    TY_SolveEquations(z1, z2, k1, k2, volumefraction, &a, &b, &c1, &c2, &d1, &d2,debug );
    checkFlag = TY_CheckSolution( z1, z2, k1, k2, volumefraction, a, b, c1, c2, d1, d2 );

    return SqTwoYukawa( qexp * 2 * radius, z1, z2, k1, k2, volumefraction, a, b, c1, c2, d1, d2 );

}