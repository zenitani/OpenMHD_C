#include <stdio.h>
#include <stdlib.h>

#include "param.hpp"
#include "u2v.hpp"
#include "set_dt.hpp"
#include "rk.hpp"
#include "fileio.hpp"
#include "limiter.hpp"
#include "glm.hpp"
// prototypes
void model( double **U, double **V, double *x, double *y, double *dx, const int ix, const int jx);
void bc( double **U, const int ix, const int jx );
void flux_solver(double **F, double **VL, double **VR,
		 const double ch, const int ix, const int jx, const int dir, const int type);
// end prototypes

using namespace param;

int main( int argc, char **argv ) {
  
  const int ix = 200 + 2;
  const int jx = 200 + 2;
  const int ijx = ix * jx;
  const int loop_max = 30000;
  const double tend  = 4.0;
  const double dtout = 0.1;
  const double cfl = 0.4;
  // Slope limiter  (0: flat, 1: minmod, 2: MC, 3: van Leer, 4: Koren)
  const int lm_type   = 1;
  // Numerical flux (0: LLF, 1: HLL, 2: HLLC, 3: HLLD)
  const int flux_type = 3;
  // Time marching  (0: TVD RK2, 1: RK2)
  const int time_type = 0;
  double t  = 0.0;
  double dt = 0.1;
  double dx;

  double x[ix], y[jx];
  
  int    i,j,k;
  int    n_output;
  double t_output;
  double ch;
  char filename[256];
  // ------------------------------------------------------------------------

  // http://www.natural-science.or.jp/article/20090911123404.php
  double **U  = new double*[var1]; // conserved variables (U)
  double **U1 = new double*[var1]; // conserved variables: medium state (U*)
  double **V  = new double*[var2]; // primitive variables (V)
  double **VR = new double*[var1];
  double **VL = new double*[var1];
  double **F  = new double*[var1];
  double **G  = new double*[var1];
  for(k=0; k<var1; k++){
    U[k]  = new double[ijx]; U1[k] = new double[ijx];
    VR[k] = new double[ijx]; VL[k] = new double[ijx];
    F[k]  = new double[ijx]; G[k]  = new double[ijx];
  }
  for(k=0; k<var2; k++)
    V[k] = new double[ijx];


  model(U,V,x,y,&dx,ix,jx);
  bc(U,ix,jx);
  set_dt(U,V,&ch,&dt,dx,cfl,ix,jx);

  t_output = -dt / 3.0;
  n_output = 0;

  if ( dt > dtout ) {
    fprintf( stderr, "error: %e > %e\n", dt, dtout);
    exit(-1);
  }
  printf("[Params]\n");
  printf(" dt: %10.3e dtout: %10.3e grids: %5d %5d\n", dt, dtout, ix, jx);
  printf(" limiter: %d  flux: %d  time-marching: %d\n", lm_type, flux_type, time_type );
  printf("== start ==\n");

  for( int k=1; k<=loop_max; k++ ){

    printf("t = %f\n", t );
    u2v(U,V,ix,jx);

    if ( t >= t_output ) {
      printf( "data output   t = %f\n", t );
      sprintf( filename, "data/field-%05d.dat", n_output );
      //printf( "%s\n", filename );
      fileio_output(filename,t,x,y,U,V,ix,jx);
      n_output++;
      t_output += dtout;
    }
    if ( t >= tend ) break;
    if ( k >= loop_max ) {
      fprintf( stderr, "max loop\n" );
      break;
    }
    set_dt(U,V,&ch,&dt,dx,cfl,ix,jx);
    glm_ss2(U,ch,dt,ix,jx);

    limiter(V[vx],VL[vx],VR[vx],ix,jx,1,lm_type);
    limiter(V[vy],VL[vy],VR[vy],ix,jx,1,lm_type);
    limiter(V[vz],VL[vz],VR[vz],ix,jx,1,lm_type);
    limiter(V[pr],VL[pr],VR[pr],ix,jx,1,lm_type);
    limiter(U[ro],VL[ro],VR[ro],ix,jx,1,lm_type);
    limiter(U[bx],VL[bx],VR[bx],ix,jx,1,lm_type);
    limiter(U[by],VL[by],VR[by],ix,jx,1,lm_type);
    limiter(U[bz],VL[bz],VR[bz],ix,jx,1,lm_type);
    limiter(U[ps],VL[ps],VR[ps],ix,jx,1,lm_type);

    for(k=0; k<var1; k++){
      for(j=0; j<jx; j++){
 	VR[k][j*ix+ix-2] = VR[k][j*ix+0];
 	VL[k][j*ix+0]    = VL[k][j*ix+ix-2];
      }
    }
    flux_solver(F,VL,VR,ch,ix,jx,1,flux_type);

    limiter(V[vx],VL[vx],VR[vx],ix,jx,2,lm_type);
    limiter(V[vy],VL[vy],VR[vy],ix,jx,2,lm_type);
    limiter(V[vz],VL[vz],VR[vz],ix,jx,2,lm_type);
    limiter(V[pr],VL[pr],VR[pr],ix,jx,2,lm_type);
    limiter(U[ro],VL[ro],VR[ro],ix,jx,2,lm_type);
    limiter(U[bx],VL[bx],VR[bx],ix,jx,2,lm_type);
    limiter(U[by],VL[by],VR[by],ix,jx,2,lm_type);
    limiter(U[bz],VL[bz],VR[bz],ix,jx,2,lm_type);
    limiter(U[ps],VL[ps],VR[ps],ix,jx,2,lm_type);

    for(k=0; k<var1; k++){
      for(i=0; i<ix; i++){
 	VL[k][i]           = VL[k][(jx-2)*ix+i];
 	VR[k][(jx-2)*ix+i] = VR[k][i];
      }
    }
    flux_solver(G,VL,VR,ch,ix,jx,2,flux_type);

    rk_tvd21(U,U1,F,G,dt,dx,ix,jx);

    bc(U1,ix,jx);
    u2v(U1,V,ix,jx);

    limiter(V[vx],VL[vx],VR[vx],ix,jx,1,lm_type);
    limiter(V[vy],VL[vy],VR[vy],ix,jx,1,lm_type);
    limiter(V[vz],VL[vz],VR[vz],ix,jx,1,lm_type);
    limiter(V[pr],VL[pr],VR[pr],ix,jx,1,lm_type);
    limiter(U1[ro],VL[ro],VR[ro],ix,jx,1,lm_type);
    limiter(U1[bx],VL[bx],VR[bx],ix,jx,1,lm_type);
    limiter(U1[by],VL[by],VR[by],ix,jx,1,lm_type);
    limiter(U1[bz],VL[bz],VR[bz],ix,jx,1,lm_type);
    limiter(U1[ps],VL[ps],VR[ps],ix,jx,1,lm_type);

    for(k=0; k<var1; k++){
      for(j=0; j<jx; j++){
 	VR[k][j*ix+ix-2] = VR[k][j*ix+0];
 	VL[k][j*ix+0]    = VL[k][j*ix+ix-2];
      }
    }
    flux_solver(F,VL,VR,ch,ix,jx,1,flux_type);

    limiter(V[vx],VL[vx],VR[vx],ix,jx,2,lm_type);
    limiter(V[vy],VL[vy],VR[vy],ix,jx,2,lm_type);
    limiter(V[vz],VL[vz],VR[vz],ix,jx,2,lm_type);
    limiter(V[pr],VL[pr],VR[pr],ix,jx,2,lm_type);
    limiter(U1[ro],VL[ro],VR[ro],ix,jx,2,lm_type);
    limiter(U1[bx],VL[bx],VR[bx],ix,jx,2,lm_type);
    limiter(U1[by],VL[by],VR[by],ix,jx,2,lm_type);
    limiter(U1[bz],VL[bz],VR[bz],ix,jx,2,lm_type);
    limiter(U1[ps],VL[ps],VR[ps],ix,jx,2,lm_type);

    for(k=0; k<var1; k++){
      for(i=0; i<ix; i++){
 	VL[k][i]           = VL[k][(jx-2)*ix+i];
 	VR[k][(jx-2)*ix+i] = VR[k][i];
      }
    }
    flux_solver(G,VL,VR,ch,ix,jx,2,flux_type);

    rk_tvd22(U,U1,F,G,dt,dx,ix,jx);

    glm_ss2(U,ch,dt,ix,jx);
    bc(U,ix,jx);
    t += dt;

  }

  return 0;
}
