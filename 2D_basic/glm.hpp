#include <math.h>
#include "param.hpp"
using namespace param;

void glm_ss2( double **U, const double ch, const double dt, const int ix, const int jx)
{
  const double cr = 0.18;
  const double f1 = exp( -0.5*dt/ch/cr );
  const int ijx = ix * jx;

  for(int i=0; i<ijx; i++)
    U[ps][i] *= f1;
  
  return;
}
