#include <stdio.h>
#include <stdlib.h>
#include "param.hpp"

using namespace param;

void set_dt( double **U, double **V, double *vmax, double *dt,
	     const double dx, const double cfl, const int ix, const int jx)
{
  int i, j, ij;  
  const double dtmin = 1.0e-7;
  double vtmp[ix*jx];
  double B2, f1, f2, vfx, vfy;
  const int is = MIN(1,ix-1), ie = MAX(1,ix-1);
  const int js = MIN(1,jx-1), je = MAX(1,jx-1);
  *vmax = -1.0;

  for(j=js; j<je; j++){
    for(i=is; i<ie; i++){

      ij = ix*j + i;
      B2 = U[bx][ij]*U[bx][ij] + U[by][ij]*U[by][ij] + U[bz][ij]*U[bz][ij];
      f1 = param::gamma * V[pr][ij];

      // fast mode in the X direction
      f2 = 4 * f1 * U[bx][ij] * U[bx][ij];
      vfx = sqrt( ( (f1+B2) + sqrt(MAX( (f1+B2)*(f1+B2)-f2, 0.0 ))) / ( 2*U[ro][ij] ));
      // fast mode in the Y direction
      f2 = 4 * f1 * U[by][ij] * U[by][ij];
      vfy = sqrt( ( (f1+B2) + sqrt(MAX( (f1+B2)*(f1+B2)-f2, 0.0 ))) / ( 2*U[ro][ij] ));

      // max speed of MHD waves
      vtmp[ij] = MAX( ABS( V[vx][ij] ) + vfx, ABS( V[vy][ij] + vfy ));
      *vmax = MAX( vtmp[ij], *vmax );
    }
  }

  *dt = cfl * dx / (*vmax);
  
  if( *dt <= dtmin ){
    fprintf( stderr, "Dt is too small.\n" );
    exit(-1);
  }

  return;
}
