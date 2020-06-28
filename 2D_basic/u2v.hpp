#include "param.hpp"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace param;

void u2v( double **U, double **V, const int ix, const int jx )
{
  const double f1 = param::gamma - 1.0;
  const int ijx = ix * jx;
  double rv2, B2;
  double prmin  = U[en][0];
  double rhomin = U[ro][0];

  for(int i=0; i<ijx; i++){
    V[vx][i] = U[mx][i] / U[ro][i];
    V[vy][i] = U[my][i] / U[ro][i];
    V[vz][i] = U[mz][i] / U[ro][i];

    rv2 = V[vx][i]*U[mx][i] + V[vy][i]*U[my][i] + V[vz][i]*U[mz][i];
    B2  = U[bx][i]*U[bx][i] + U[by][i]*U[by][i] + U[bz][i]*U[bz][i];
    V[pr][i] = f1 * ( U[en][i] - 0.5 * ( rv2 + B2 ));

    prmin  = MIN( V[pr][i], prmin );
    rhomin = MIN( U[ro][i], rhomin );
  }
  
  if( prmin <= 0.0 ){
    fprintf( stderr, "Negative pressure at somewhere: %e\n", prmin );
    exit(-1);
  }
  if( rhomin <= 0.0 ){
    fprintf( stderr, "Negative density at somewhere: %e\n", rhomin );
    exit(-1);
  }

  return;
}
