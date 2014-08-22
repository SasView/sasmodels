real J1(real x);
real J1(real x)
{
  const real ax = fabs(x);
  if (ax < REAL(8.0)) {
    const real y = x*x;
    const real ans1 = x*(REAL(72362614232.0)
              + y*(REAL(-7895059235.0)
              + y*(REAL(242396853.1)
              + y*(REAL(-2972611.439)
              + y*(REAL(15704.48260)
              + y*(REAL(-30.16036606)))))));
    const real ans2 = REAL(144725228442.0)
              + y*(REAL(2300535178.0)
              + y*(REAL(18583304.74)
              + y*(REAL(99447.43394)
              + y*(REAL(376.9991397)
              + y))));
    return ans1/ans2;
  } else {
    const real y = REAL(64.0)/(ax*ax);
    const real xx = ax - REAL(2.356194491);
    const real ans1 = REAL(1.0)
              + y*(REAL(0.183105e-2)
              + y*(REAL(-0.3516396496e-4)
              + y*(REAL(0.2457520174e-5)
              + y*REAL(-0.240337019e-6))));
    const real ans2 = REAL(0.04687499995)
              + y*(REAL(-0.2002690873e-3)
              + y*(REAL(0.8449199096e-5)
              + y*(REAL(-0.88228987e-6)
              + y*REAL(0.105787412e-6))));
    real sn,cn;
    SINCOS(xx, sn, cn);
    const real ans = sqrt(REAL(0.636619772)/ax) * (cn*ans1 - (REAL(8.0)/ax)*sn*ans2);
    return (x < REAL(0.0)) ? -ans : ans;
  }
}
