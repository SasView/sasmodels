double J1(double x);
double J1(double x)
{
  const double ax = fabs(x);
  if (ax < 8.0) {
    const double y = x*x;
    const double ans1 = x*(72362614232.0
              + y*(-7895059235.0
              + y*(242396853.1
              + y*(-2972611.439
              + y*(15704.48260
              + y*(-30.16036606))))));
    const double ans2 = 144725228442.0
              + y*(2300535178.0
              + y*(18583304.74
              + y*(99447.43394
              + y*(376.9991397
              + y))));
    return ans1/ans2;
  } else {
    const double y = 64.0/(ax*ax);
    const double xx = ax - 2.356194491;
    const double ans1 = 1.0
              + y*(0.183105e-2
              + y*(-0.3516396496e-4
              + y*(0.2457520174e-5
              + y*-0.240337019e-6)));
    const double ans2 = 0.04687499995
              + y*(-0.2002690873e-3
              + y*(0.8449199096e-5
              + y*(-0.88228987e-6
              + y*0.105787412e-6)));
    double sn,cn;
    SINCOS(xx, sn, cn);
    const double ans = sqrt(0.636619772/ax) * (cn*ans1 - (8.0/ax)*sn*ans2);
    return (x < 0.0) ? -ans : ans;
  }
}
