#include "param.hpp"
#include <stdio.h>
#include <stdlib.h>

using namespace param;

void fileio_output( char filename[], double t, double x[], double y[],
		    double **U, double **V, const int ix, const int jx )
{
  FILE *fp;
  const int ijx = ix*jx;

  if(( fp=fopen(filename,"wb") )==NULL ){
    fprintf( stderr, "cannot open file: %s\n", filename );
    exit(-1);
  }
  fwrite(&t,sizeof(double),1,fp);
  fwrite(&ix,sizeof(int),1,fp);
  fwrite(&jx,sizeof(int),1,fp);
  fwrite(x,sizeof(double),ix,fp);
  fwrite(y,sizeof(double),jx,fp);
  fwrite(U[mx],sizeof(double),ijx,fp);
  fwrite(U[my],sizeof(double),ijx,fp);
  fwrite(U[mz],sizeof(double),ijx,fp);
  fwrite(U[en],sizeof(double),ijx,fp);
  fwrite(U[ro],sizeof(double),ijx,fp);
  fwrite(U[bx],sizeof(double),ijx,fp);
  fwrite(U[by],sizeof(double),ijx,fp);
  fwrite(U[bz],sizeof(double),ijx,fp);
  fwrite(U[ps],sizeof(double),ijx,fp);
  fwrite(V[vx],sizeof(double),ijx,fp);
  fwrite(V[vy],sizeof(double),ijx,fp);
  fwrite(V[vz],sizeof(double),ijx,fp);
  fwrite(V[pr],sizeof(double),ijx,fp);
  fclose(fp);
  return;
}
