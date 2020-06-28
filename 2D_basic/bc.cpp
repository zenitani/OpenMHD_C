#include "param.hpp"

using namespace param;

void bc( double **U, const int ix, const int jx )
{
  int i,j,k;
  for(k=0; k<var1; k++){
    for(j=0; j<jx; j++){
      U[k][INDEX(j,ix-1)] = U[k][INDEX(j,   1)];
      U[k][INDEX(j,   0)] = U[k][INDEX(j,ix-2)];
    }
  }  
  for(k=0; k<var1; k++){
    for(i=0; i<ix; i++){
      U[k][INDEX(jx-1,i)] = U[k][INDEX(   1,i)];
      U[k][INDEX(   0,i)] = U[k][INDEX(jx-2,i)];
    }
  }

  return;
}
