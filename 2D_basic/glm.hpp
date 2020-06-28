#include <math.h>
#include "param.hpp"
using namespace param;

void glm_ss( double **U, const double ch, const double dt, const int ix, const int jx)
{
  const double cr = 0.18;
  const double f1 = exp( -dt/ch/cr );
  const int ijx = ix * jx;

  for(int i=0; i<ijx; i++)
    U[ps][i] *= f1;
  
  return;
}

void flux_glm(double **F, double **VL, double **VR,
	      const int ch, const int ix, const int jx, const int dir)
{
  const double ch2 = ch*ch;
  int i, j, ij, is=0, ie=ix-1, js=0, je=jx-1;
  int bn = bx + (dir-1) % 3;

  switch(dir){
  case 1:
    js = 1;
    break;
  case 2:
    is = 1;
    break;
  case 3:
    break;
  }

  for(j=js;j<je;j++){
    for(i=is;i<ie;i++){
      ij = j*ix + i;
      F[bn][ij] = 0.5*( VL[ps][ij] + VR[ps][ij] - ch*( VR[bn][ij]-VL[bn][ij] ) );
      F[ps][ij] = 0.5*( ch2*( VL[bn][ij]+VR[bn][ij] ) - ch*( VR[ps][ij]-VL[ps][ij] ) );
    }
  }

  return;
}

