#include "param.hpp"

void rk_tvd21(double **wk, double **wk1, double **wF, double **wG,
	      const double dt, const double dx, const int ix, const int jx)
{
  const double dtx = dt/dx;
  int i,j,k,ij;

  for(k=0; k<var1; k++){
    for(j=1; j<jx-1; j++){
      for(i=1; i<ix-1; i++){
	ij = j*ix + i;
	wk1[k][ij] = wk[k][ij] + 
	  dtx * ( wF[k][INDEX(j,i-1)] - wF[k][ij] + wG[k][INDEX(j-1,i)] - wG[k][ij] );
      }
    }
  }  

  return;
}
void rk_tvd22(double **wk, double **wk1, double **wF, double **wG,
	      const double dt, const double dx, const int ix, const int jx)
{
  const double dtx = dt/dx;
  int i,j,k,ij;

  for(k=0; k<var1; k++){
    for(j=1; j<jx-1; j++){
      for(i=1; i<ix-1; i++){
	ij = j*ix + i;
	wk[k][ij] = 0.5*( wk[k][ij] + wk1[k][ij] +
			  dtx * ( wF[k][INDEX(j,i-1)] - wF[k][ij] + wG[k][INDEX(j-1,i)] - wG[k][ij] ) );
      }
    }
  }  

  return;
}
