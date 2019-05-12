#include <math.h> // log(), sqrt()

#include "utils.h" // rand01()
#include "md_utils.h"


double gauss(double stdev, double mean)
// gaussian distribution
// see 579 page [Daan Frenkel]:
//  generate value, that gaussian distributed with mean = mean, and standart deviation = stdev
{
  double v1, v2, res;
  double r = 2.0;

  while (r > 1.0)
  {
    v1 = 2.0 * rand01() - 1.0;
    v2 = 2.0 * rand01() - 1.0;
    r = v1 * v1 + v2 * v2;
    //printf("v1=%f;    v2=%f;    r=%f\n", v1, v2, r);
  }
  res = v1 * sqrt(-2.0*log(r)/r);
  res = mean + stdev * res;
  return res;
}
// end gauss function

