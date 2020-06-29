#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "param.hpp"
inline double sign(double a, double b) { return( (b < 0.0) ? (-a) : a ); }

using namespace param;

void limiter( double *wk, double *wL, double *wR,
	      const int ix, const int jx, const int dir, const int type)
{
  int i,j,ij, is=0,ie=0,js=0,je=0;
  double gA, gB, grad;
  const int ijx = ix*jx;

  for(i=0; i<ijx; i++){
    wL[i] = 0.0;
    wR[i] = 0.0;
  }

  switch( dir ){
  case 1:
    js = MIN(1,jx-1);
    je = MAX(1,jx-1);
    break;
  case 2:
    is = MIN(1,ix-1);
    ie = MAX(1,ix-1);
    break;
  }

  switch( type ){
  case 1:

    switch( dir ){
    case 1:
      for(j=js; j<je; j++){
	wL[INDEX(j,0)] = wk[INDEX(j,0)];
	for(i=1; i<ix-1; i++){
	  ij = j*ix+i;
	  gA = wk[ij]   - wk[ij-1];
	  gB = wk[ij+1] - wk[ij];
	  grad = (sign(0.25,gA)+sign(0.25,gB))*MIN(fabs(gA),fabs(gB));
	  wR[ij-1] = wk[ij] - grad;
	  wL[ij]   = wk[ij] + grad;
	}
	wR[INDEX(j,ix-2)] = wk[INDEX(j,ix-1)];
      }
      break;

    case 2:
      for(i=is; i<ie; i++){
	wL[i] = wk[i];
      }
      for(j=1; j<jx-1; j++){
	for(i=is; i<ie; i++){
	  gA = wk[INDEX(j,i)]   - wk[INDEX(j-1,i)];
	  gB = wk[INDEX(j+1,i)] - wk[INDEX(j,i)];
	  grad = (sign(0.25,gA)+sign(0.25,gB))*MIN(fabs(gA),fabs(gB));
	  wR[INDEX(j-1,i)] = wk[INDEX(j,i)] - grad;
	  wL[INDEX(j,i)]   = wk[INDEX(j,i)] + grad;
	}
      }
      for(i=is; i<ie; i++){
	wR[INDEX(jx-2,i)] = wk[INDEX(jx-2,i)];
      }
      break;
    }
    break;

  default:
    fprintf( stderr, "unknown limiter\n" );
    exit(-1);
    break;
  }

  return;
}
