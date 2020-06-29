#include <math.h>
#include "param.hpp"
using namespace param;

void flux_solver(double **F, double **VL, double **VR,
		 const int ix, const int jx, const int dir, const int type)
{
  const double EPS = 1e-8;
  const double r1 = 1.0 / ( param::gamma - 1.0 );
  //const double r2 = param::gamma / ( param::gamma - 1.0 );
  const int ijx = ix*jx;
  int i,j,k,ij;
  int is=0, ie=ix-1, js=0, je=jx-1;
  int vn, vt1, vt2, mn, mt1, mt2, bn, bt1, bt2;

  for(k=0;k<var1;k++)
    for(i=0;i<ijx;i++)
      F[k][i] = 0.0;

  // directions
  vn  = vx + (dir-1) % 3;
  vt1 = vx + (dir  ) % 3;
  vt2 = vx + (dir+1) % 3;
  bn  = bx + (dir-1) % 3;
  bt1 = bx + (dir  ) % 3;
  bt2 = bx + (dir+1) % 3;
  mn  = mx + (dir-1) % 3;
  mt1 = mx + (dir  ) % 3;
  mt2 = mx + (dir+1) % 3;

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

  switch(type){
  case 3:

    double rol,vnl,vtl,vul,btl,bul,prl,ror,vnr,vtr,vur,btr,bur,prr,bnc;
    for(j=js;j<je;j++){
      for(i=is;i<ie;i++){
	ij = j*ix + i;
	rol = VL[ro][ij];
	vnl = VL[vn][ij];
	vtl = VL[vt1][ij];
	vul = VL[vt2][ij];
	btl = VL[bt1][ij];
	bul = VL[bt2][ij];
	prl = VL[pr][ij];
	bnc = 0.5*(VL[bn][ij]+VR[bn][ij]);
	ror = VR[ro][ij];
	vnr = VR[vn][ij];
	vtr = VR[vt1][ij];
	vur = VR[vt2][ij];
	btr = VR[bt1][ij];
	bur = VR[bt2][ij];
	prr = VR[pr][ij];

	/* Bn at the interface */
	double bnc2=bnc*bnc;
	int sgn=(bnc > 0)?(1):(-1);
	/* Variables at the left-face */
	double roli=1.0/rol;
	double pml=0.5*(btl*btl+bul*bul);
	double ptl=prl+pml;
	double enl=r1*prl+pml+0.5*rol*(vnl*vnl+vtl*vtl+vul*vul);
	double vbl=vtl*btl+vul*bul;
	/* Variables at the right-face */
	double rori=1.0/ror;
	double pmr=0.5*(btr*btr+bur*bur);
	double ptr=prr+pmr;
	double enr=r1*prr+pmr+0.5*ror*(vnr*vnr+vtr*vtr+vur*vur);
	double vbr=vtr*btr+vur*bur;
	/* Maximum/minimum wave speeds */
	double cl2=param::gamma*prl*roli;
	double cr2=param::gamma*prr*rori;
	double cal2=bnc2*roli;
	double car2=bnc2*rori;
	double cbl2=cl2+cal2+2.0*pml*roli;
	double cbr2=cr2+car2+2.0*pmr*rori;
	double cfl2=0.5*(cbl2+sqrt(fabs(cbl2*cbl2-4.0*cl2*cal2)));
	double cfr2=0.5*(cbr2+sqrt(fabs(cbr2*cbr2-4.0*cr2*car2)));
	double cfl=sqrt(cfl2);
	double cfr=sqrt(cfr2);
	double sl=MIN(0.0,MIN(vnl,vnr)-MAX(cfl,cfr));
	double sr=MAX(0.0,MAX(vnl,vnr)+MAX(cfl,cfr));
	/* HLL average of the normal velocity and the total pressure */
	double slvl=sl-vnl;
	double srvr=sr-vnr;
	double rslvl=rol*slvl;
	double rsrvr=ror*srvr;
	double drsvi=1.0/(rsrvr-rslvl);
	double vnc=(rsrvr*vnr-rslvl*vnl-ptr+ptl)*drsvi;
	double ptc=(rsrvr*ptl-rslvl*ptr+rsrvr*rslvl*(vnr-vnl))*drsvi;
	/* Variables of the outer sides in the Riemann fan */
	double slvc=sl-vnc;
	double srvc=sr-vnc;
	double ro2l=rslvl/slvc;
	double ro2r=rsrvr/srvc;
	double rhdl=rslvl*slvc-bnc2;
	double rhdr=rsrvr*srvc-bnc2;
	double vt2l,vu2l,bt2l,bu2l;
	double vt2r,vu2r,bt2r,bu2r;
	if (fabs(rhdl) > EPS){
	  double rhdli=1.0/rhdl;
	  double rhnvl=(vnl-vnc)*bnc;
	  double rhnbl=rslvl*slvl-bnc2;
	  vt2l=vtl+rhnvl*rhdli*btl;
	  vu2l=vul+rhnvl*rhdli*bul;
	  bt2l=rhnbl*rhdli*btl;
	  bu2l=rhnbl*rhdli*bul;
	} else{
	  vt2l=vtl;
	  vu2l=vul;
	  bt2l=btl;
	  bu2l=bul;
	}
	if (fabs(rhdr) > EPS){
	  double rhdri=1.0/rhdr;
	  double rhnvr=(vnr-vnc)*bnc;
	  double rhnbr=rsrvr*srvr-bnc2;
	  vt2r=vtr+rhnvr*rhdri*btr;
	  vu2r=vur+rhnvr*rhdri*bur;
	  bt2r=rhnbr*rhdri*btr;
	  bu2r=rhnbr*rhdri*bur;
	} else{
	  vt2r=vtr;
	  vu2r=vur;
	  bt2r=btr;
	  bu2r=bur;
	}
	double vb2l=vt2l*bt2l+vu2l*bu2l;
	double vb2r=vt2r*bt2r+vu2r*bu2r;
	double en2l=(slvl*enl-ptl*vnl+ptc*vnc+bnc*(vbl-vb2l))/slvc;
	double en2r=(srvr*enr-ptr*vnr+ptc*vnc+bnc*(vbr-vb2r))/srvc;
	/* Variables of the inner sides in the Riemann fan */
	double rro2l=sqrt(ro2l);
	double rro2r=sqrt(ro2r);
	double rro2i=1.0/(rro2r+rro2l);
	double vt3m=(rro2r*vt2r+rro2l*vt2l+(bt2r-bt2l)*sgn)*rro2i;
	double vu3m=(rro2r*vu2r+rro2l*vu2l+(bu2r-bu2l)*sgn)*rro2i;
	double bt3m=(rro2l*bt2r+rro2r*bt2l+rro2r*rro2l*(vt2r-vt2l)*sgn)*rro2i;
	double bu3m=(rro2l*bu2r+rro2r*bu2l+rro2r*rro2l*(vu2r-vu2l)*sgn)*rro2i;
	double vb3m=vt3m*bt3m+vu3m*bu3m;
	double en3l=en2l-rro2l*(vb2l-vb3m)*sgn;
	double en3r=en2r+rro2r*(vb2r-vb3m)*sgn;
	/* Variables at the interface */
	double rou,vtu,vuu,btu,buu,enu;
	if (vnc-fabs(bnc)/rro2l > 0){
	  rou=ro2l;
	  vtu=vt2l;
	  vuu=vu2l;
	  btu=bt2l;
	  buu=bu2l;
	  enu=en2l;
	} else{
	  if (vnc >= 0){
	    rou=ro2l;
	    vtu=vt3m;
	    vuu=vu3m;
	    btu=bt3m;
	    buu=bu3m;
	    enu=en3l;
	  } else{
	    if (vnc+fabs(bnc)/rro2r >= 0){
	      rou=ro2r;
	      vtu=vt3m;
	      vuu=vu3m;
	      btu=bt3m;
	      buu=bu3m;
	      enu=en3r;
	    } else{
	      rou=ro2r;
	      vtu=vt2r;
	      vuu=vu2r;
	      btu=bt2r;
	      buu=bu2r;
	      enu=en2r;
	    }
	  }
	}
	/* HLLD fluxes */
	F[ro][ij] = rou*vnc;
	F[mn][ij] = rou*vnc*vnc+ptc-bnc2*0.5;
	F[mt1][ij]= rou*vtu*vnc-bnc*btu;
	F[mt2][ij]= rou*vuu*vnc-bnc*buu;
	F[bt1][ij]= btu*vnc-bnc*vtu;
	F[bt2][ij] = buu*vnc-bnc*vuu;
	F[en][ij] = (enu+ptc)*vnc-bnc*(vtu*btu+vuu*buu);

      }
    }

    break;
  }
}
