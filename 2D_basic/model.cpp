#include "param.hpp"
#include <math.h>
#include <stdio.h>
using namespace param;

void model( double **U, double **V, double *x, double *y, double *dx, const int ix, const int jx )
{
  const double domain_x[2] = {0.0, 2*pi};
  const double domain_y[1] = {0.0};
  const double f1 = 1.0 / ( param::gamma - 1 );
  double v2, B2;
  int i, j, ij;

  *dx = ( domain_x[1] - domain_x[0] ) / double( ix-2 );
  x[0] = domain_x[0] - (*dx)/2;
  for(i=1; i<ix; i++)
    x[i] = x[i-1] + (*dx);
  y[0] = domain_y[0] - (*dx)/2;
  for(j=1; j<jx; j++)
    y[j] = y[j-1] + (*dx);

  for(j=0; j<jx; j++){
    for(i=0; i<ix; i++){

      ij = ix*j + i;
      U[ro][ij] = param::gamma * param::gamma;

      V[pr][ij] = param::gamma;
      V[vx][ij] = -sin(y[j]);
      V[vy][ij] = sin(x[i]);
      V[vz][ij] = 0.0;

      U[bx][ij] = -sin(y[j]);
      U[by][ij] = sin(2*x[i]);
      U[bz][ij] = 0.0;
      U[ps][ij] = 0.0;

      v2 = V[vx][ij]*V[vx][ij] + V[vy][ij]*V[vy][ij] + V[vz][ij]*V[vz][ij];
      B2 = U[bx][ij]*U[bx][ij] + U[by][ij]*U[by][ij] + U[bz][ij]*U[bz][ij];
      U[mx][ij] = U[ro][ij] * V[vx][ij];
      U[my][ij] = U[ro][ij] * V[vy][ij];
      U[mz][ij] = U[ro][ij] * V[vz][ij];
      U[en][ij] = 0.5 * ( U[ro][ij]*v2 + B2 ) + f1*V[pr][ij];
    }
  }
  
  return;
}
