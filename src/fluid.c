#include "run_param.h"
#include "fluid.h"

void calc_flux_jacobian_plus(REAL *Wp, REAL *L, REAL *R, REAL *lambda)
{
#define LL(a,b) L[(a)+3*(b)]  
#define RR(a,b) R[(a)+3*(b)]
#define WP(a,b) Wp[(a)+3*(b)]
  for(int ii=0;ii<3;ii++) {
    for(int jj=0;jj<3;jj++) {
      WP(ii,jj) = 0.0; 
      for(int m=0;m<3;m++) WP(ii,jj) += RR(ii,m)*fmax(lambda[m],0.0)*LL(m,jj);  
    }
  }
#undef LL
#undef RR
#undef WP
}

void calc_flux_jacobian_minus(REAL *Wm, REAL *L, REAL *R, REAL *lambda)
{
#define LL(a,b) L[(a)+3*(b)]  
#define RR(a,b) R[(a)+3*(b)]
#define WM(a,b) Wm[(a)+3*(b)]
  for(int ii=0;ii<3;ii++) {
    for(int jj=0;jj<3;jj++) {
      WM(ii,jj) = 0.0; 
      for(int m=0;m<3;m++) WM(ii,jj) += RR(ii,m)*fmin(lambda[m],0.0)*LL(m,jj);
    }
  }
#undef LL
#undef RR
#undef WM
}

void calc_flux_jacobian(REAL *Wp, REAL *Wm, REAL *L, REAL *R, REAL *lambda)
{
#define LL(a,b) L[(a)+3*(b)]  
#define RR(a,b) R[(a)+3*(b)]
#define WP(a,b) Wp[(a)+3*(b)]
#define WM(a,b) Wm[(a)+3*(b)]
  for(int ii=0;ii<3;ii++) {
    for(int jj=0;jj<3;jj++) {
      WP(ii,jj) = 0.0; 
      WM(ii,jj) = 0.0;
      for(int m=0;m<3;m++) WP(ii,jj) += RR(ii,m)*fmax(lambda[m],0.0)*LL(m,jj);  
      for(int m=0;m<3;m++) WM(ii,jj) += RR(ii,m)*fmin(lambda[m],0.0)*LL(m,jj);
    }
  }
#undef LL
#undef RR
#undef WP
#undef WM
}

void set_eigen_value(REAL *lambda, REAL dens, REAL velx, REAL pres)
{
  REAL cs = sqrt(GAMMA*pres/dens);
  lambda[0] = velx - cs;
  lambda[1] = velx;
  lambda[2] = velx + cs;
}

void set_LR_matrix(REAL *L, REAL *R, REAL dens, REAL velx, REAL pres)
{
  REAL cs = sqrt(GAMMA*pres/dens);
  REAL enth = 0.5*SQR(velx) + (GAMMA*pres)/((GAMMA-1.0)*dens);

#define LL(a,b) L[(a)+3*(b)]  
#define RR(a,b) R[(a)+3*(b)]
  RR(0,0) = 1.0;
  RR(0,1) = 1.0;
  RR(0,2) = 1.0;

  RR(1,0) = velx - cs;
  RR(1,1) = velx;
  RR(1,2) = velx + cs;

  RR(2,0) = enth - velx*cs;
  RR(2,1) = 0.5*SQR(velx);
  RR(2,2) = enth + velx*cs;

  REAL b1 = 0.5*(GAMMA-1.0)*SQR(velx/cs);
  REAL b2 = (GAMMA-1.0)/SQR(cs);

  LL(0,0) = 0.5*(b1 + velx/cs);
  LL(0,1) = -0.5*(1.0/cs + b2*velx);
  LL(0,2) = 0.5*b2;
  
  LL(1,0) = 1.0-b1;
  LL(1,1) = b2*velx;
  LL(1,2) = -b2;
  
  LL(2,0) = 0.5*(b1 - velx/cs);
  LL(2,1) = 0.5*(1.0/cs - b2*velx);
  LL(2,2) = 0.5*b2;
#undef LL
#undef RR
}

void convert_to_wspace(struct fluid_1d *mesh, struct run_param *this_run)
{
  static REAL L[9], R[9];
  static REAL w[3];

  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL dens = mesh[ix].dens;
    REAL velx = mesh[ix].momx/dens;
    REAL pres = (GAMMA-1.0)*(mesh[ix].eneg - 0.5*dens*SQR(velx));
    set_LR_matrix(L, R, dens, velx, pres);

#define LL(a,b) L[(a)+3*(b)]  
    for(int m=0;m<3;m++) {
      w[m] = 
	LL(m,0)*mesh[ix].dens + LL(m,1)*mesh[ix].momx + LL(m,2)*mesh[ix].eneg;
    }
#undef LL
    mesh[ix].dens = w[0];
    mesh[ix].momx = w[1];
    mesh[ix].eneg = w[2];
  }
}

int high_mach_detector(struct fluid_1d *mesh) 
{
  REAL eneg_th = (mesh->eneg - 0.5*SQR(mesh->momx)/mesh->dens);

  if(eneg_th < HIGH_MACH_THRESHOLD*mesh->eneg) {
    return 1;
  }else{
    return 0;
  }
}

int shock_detector(struct fluid_1d *mesh, int ix, struct run_param *this_run)
{
  int shock=0;
  int ixm = ix-1; if(ixm<0) ixm=0;
  int ixp = ix+1; if(ixp>this_run->nmesh-1) ixp=this_run->nmesh-1;
  
  REAL velx_ip = mesh[ixp].momx/mesh[ixp].dens;
  REAL velx_im = mesh[ixm].momx/mesh[ixm].dens;
  REAL pres_ip = 
    (GAMMA-1.0)*(mesh[ixp].eneg-0.5*SQR(mesh[ixp].momx)/mesh[ixp].dens);
  REAL pres_im = 
    (GAMMA-1.0)*(mesh[ixm].eneg-0.5*SQR(mesh[ixm].momx)/mesh[ixm].dens);
  REAL pres_i  = 
    (GAMMA-1.0)*(mesh[ix ].eneg-0.5*SQR(mesh[ix ].momx)/mesh[ix ].dens);


  if(velx_ip - velx_im < 0.0) shock = 1;
  if(fabsf(pres_ip-pres_im) > SHOCK_PRESS_THRESHOLD*pres_i) shock = 1;

  return shock;
}

void calc_flux_FVS_MP5SL(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
			 struct run_param *this_run)
{
  static REAL L[9], R[9], Wp[9], Wm[9];
  static REAL wim2[3],wim1[3],wi[3],wip1[3],wip2[3];
  static REAL wflux_iph[3],wflux_imh[3];
  static REAL lambda_L[3],del_L[3];
  static REAL lambda_R[3],del_R[3];
  static REAL lambda[3], del[3];

  REAL (*Ep)[3], (*Em)[3];

  Ep = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  Em = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixm2 = ix-2; if(ixm2<0) ixm2 = 0;
    int ixm1 = ix-1; if(ixm1<0) ixm1 = 0;
    int ixp1 = ix+1; if(ixp1>this_run->nmesh-1) ixp1 = this_run->nmesh-1;
    int ixp2 = ix+2; if(ixp2>this_run->nmesh-1) ixp2 = this_run->nmesh-1;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    REAL dens_im2 = mesh[ixm2].dens;
    REAL dens_im1 = mesh[ixm1].dens;
    REAL dens_i   = mesh[ix  ].dens;
    REAL dens_ip1 = mesh[ixp1].dens;
    REAL dens_ip2 = mesh[ixp2].dens;

    REAL velx_im2 = mesh[ixm2].momx/mesh[ixm2].dens;
    REAL velx_im1 = mesh[ixm1].momx/mesh[ixm1].dens;
    REAL velx_i   = mesh[ix].momx/mesh[ix].dens;
    REAL velx_ip1 = mesh[ixp1].momx/mesh[ixp1].dens;
    REAL velx_ip2   = mesh[ixp2].momx/mesh[ixp2].dens;

    REAL pres_im2 = (GAMMA-1)*(mesh[ixm2].eneg - 0.5*dens_im2*SQR(velx_im2));
    REAL pres_im1 = (GAMMA-1)*(mesh[ixm1].eneg - 0.5*dens_im1*SQR(velx_im1));
    REAL pres_i   = (GAMMA-1)*(mesh[ix  ].eneg - 0.5*dens_i  *SQR(velx_i  ));
    REAL pres_ip1 = (GAMMA-1)*(mesh[ixp1].eneg - 0.5*dens_ip1*SQR(velx_ip1));
    REAL pres_ip2 = (GAMMA-1)*(mesh[ixp2].eneg - 0.5*dens_ip2*SQR(velx_ip2));

    REAL dens_L = mp5_L(dens_im2, dens_im1, dens_i, dens_ip1, dens_ip2);
    REAL velx_L = mp5_L(velx_im2, velx_im1, velx_i, velx_ip1, velx_ip2);
    REAL pres_L = mp5_L(pres_im2, pres_im1, pres_i, pres_ip1, pres_ip2);

    REAL dens_R = mp5_R(dens_im2, dens_im1, dens_i, dens_ip1, dens_ip2);
    REAL velx_R = mp5_R(velx_im2, velx_im1, velx_i, velx_ip1, velx_ip2);
    REAL pres_R = mp5_R(pres_im2, pres_im1, pres_i, pres_ip1, pres_ip2);

    set_LR_matrix(L, R, dens_i, velx_i, pres_i);

#define LL(a,b) L[(a)+3*(b)]  
    for(int m=0;m<3;m++) {
      wim2[m] = LL(m,0)*mesh[ixm2].dens 
	      + LL(m,1)*mesh[ixm2].momx 
	      + LL(m,2)*mesh[ixm2].eneg;
      wim1[m] = LL(m,0)*mesh[ixm1].dens 
	      + LL(m,1)*mesh[ixm1].momx 
	      + LL(m,2)*mesh[ixm1].eneg;
      wi[m]   = LL(m,0)*mesh[ix  ].dens 
	      + LL(m,1)*mesh[ix  ].momx 
	      + LL(m,2)*mesh[ix  ].eneg;
      wip1[m] = LL(m,0)*mesh[ixp1].dens 
	      + LL(m,1)*mesh[ixp1].momx 
	      + LL(m,2)*mesh[ixp1].eneg;
      wip2[m] = LL(m,0)*mesh[ixp2].dens 
	      + LL(m,1)*mesh[ixp2].momx 
	      + LL(m,2)*mesh[ixp2].eneg;
    }
#undef LL
#if 1
    set_eigen_value(lambda, dens_i, velx_i, pres_i);
    for(int m=0;m<3;m++) {
      del[m] = fabs(lambda[m])*this_run->dtime/this_run->delta_x;
    }
    for(int m=0;m<3;m++) {
      wflux_iph[m] = csl5(wim2[m], wim1[m], wi[m], wip1[m], wip2[m], 
			  this_run->delta_x, del[m]);
      wflux_iph[m] = 
	mp5_bound_L(wim2[m], wim1[m], wi[m], wip1[m], wip2[m], wflux_iph[m]);
      wflux_iph[m] *= fmax(lambda[m],0.0);

      wflux_imh[m] = csl5(wip2[m], wip1[m], wi[m], wim1[m], wim2[m],
			  this_run->delta_x, del[m]);
      wflux_imh[m] = 
	mp5_bound_R(wim2[m], wim1[m], wi[m], wip1[m], wip2[m], wflux_imh[m]);
      wflux_imh[m] *= fmin(lambda[m],0.0);
    }
#else
    set_eigen_value(lambda_L, dens_L, velx_L, pres_L);
    set_eigen_value(lambda_R, dens_R, velx_R, pres_R);
    for(int m=0;m<3;m++) {
      del_L[m] = fabs(lambda_L[m])*this_run->dtime/this_run->delta_x;
      del_R[m] = fabs(lambda_R[m])*this_run->dtime/this_run->delta_x;
    }
    for(int m=0;m<3;m++) {
      wflux_iph[m] = csl5(wim2[m], wim1[m], wi[m], wip1[m], wip2[m], 
			  this_run->delta_x, del_L[m]);
      wflux_iph[m] = 
	mp5_bound_L(wim2[m], wim1[m], wi[m], wip1[m], wip2[m], wflux_iph[m]);
      wflux_iph[m] *= fmax(lambda_L[m],0.0);

      wflux_imh[m] = csl5(wip2[m], wip1[m], wi[m], wim1[m], wim2[m],
			  this_run->delta_x, del_R[m]);
      wflux_imh[m] = 
	mp5_bound_R(wim2[m], wim1[m], wi[m], wip1[m], wip2[m], wflux_imh[m]);
      wflux_imh[m] *= fmin(lambda_R[m],0.0);
    }
#endif

#define RR(a,b) R[(a)+3*(b)]
    for(int m=0;m<3;m++) {
      Ep[ix][m] = 
	RR(m,0)*wflux_iph[0] + RR(m,1)*wflux_iph[1] + RR(m,2)*wflux_iph[2];
      Em[ix][m] = 
	RR(m,0)*wflux_imh[0] + RR(m,1)*wflux_imh[1] + RR(m,2)*wflux_imh[2];
    }
#undef RR
  }

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixp = ix+1;
    flux_iph[ix].dens_flux = Ep[ix][0] + Em[ixp][0];
    flux_iph[ix].momx_flux = Ep[ix][1] + Em[ixp][1];
    flux_iph[ix].eneg_flux = Ep[ix][2] + Em[ixp][2];
  }

  free(Ep);
  free(Em);
}


void calc_flux_FVS_MP5_2(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
			 struct run_param *this_run)
{
  static REAL L[9], R[9], Wp[9], Wm[9];
  static REAL wim2[3],wim1[3],wi[3],wip1[3],wip2[3];
  static REAL wiph[3],wimh[3];
  static REAL lambda[3];

  REAL (*Ep)[3], (*Em)[3];

  Ep = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  Em = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixm2 = ix-2; if(ixm2<0) ixm2 = 0;
    int ixm1 = ix-1; if(ixm1<0) ixm1 = 0;
    int ixp1 = ix+1; if(ixp1>this_run->nmesh-1) ixp1 = this_run->nmesh-1;
    int ixp2 = ix+2; if(ixp2>this_run->nmesh-1) ixp2 = this_run->nmesh-1;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);
    
    REAL dens_im2 = mesh[ixm2].dens;
    REAL dens_im1 = mesh[ixm1].dens;
    REAL dens_i   = mesh[ix].dens;
    REAL dens_ip1 = mesh[ixp1].dens;
    REAL dens_ip2 = mesh[ixp2].dens;

    REAL velx_im2 = mesh[ixm2].momx/mesh[ixm2].dens;
    REAL velx_im1 = mesh[ixm1].momx/mesh[ixm1].dens;
    REAL velx_i   = mesh[ix].momx/mesh[ix].dens;
    REAL velx_ip1 = mesh[ixp1].momx/mesh[ixp1].dens;
    REAL velx_ip2 = mesh[ixp2].momx/mesh[ixp2].dens;

    REAL pres_im2 = (GAMMA-1)*(mesh[ixm2].eneg - 0.5*dens_im2*SQR(velx_im2));
    REAL pres_im1 = (GAMMA-1)*(mesh[ixm1].eneg - 0.5*dens_im1*SQR(velx_im1));
    REAL pres_i   = (GAMMA-1)*(mesh[ix  ].eneg - 0.5*dens_i*SQR(velx_i));
    REAL pres_ip1 = (GAMMA-1)*(mesh[ixp1].eneg - 0.5*dens_ip1*SQR(velx_ip1));
    REAL pres_ip2 = (GAMMA-1)*(mesh[ixp2].eneg - 0.5*dens_ip2*SQR(velx_ip2));

    set_LR_matrix(L, R, dens_i, velx_i, pres_i);

    REAL dens_L = mp5_L(dens_im2, dens_im1, dens_i, dens_ip1, dens_ip2);
    REAL dens_R = mp5_R(dens_im2, dens_im1, dens_i, dens_ip1, dens_ip2);

    REAL velx_L = mp5_L(velx_im2, velx_im1, velx_i, velx_ip1, velx_ip2);
    REAL velx_R = mp5_R(velx_im2, velx_im1, velx_i, velx_ip1, velx_ip2);

    REAL pres_L = mp5_L(pres_im2, pres_im1, pres_i, pres_ip1, pres_ip2);
    REAL pres_R = mp5_R(pres_im2, pres_im1, pres_i, pres_ip1, pres_ip2);

    REAL momx_L = dens_L*velx_L;
    REAL momx_R = dens_R*velx_R;
    REAL eneg_L = 0.5*dens_L*SQR(velx_L) + pres_L/(GAMMA-1.0);
    REAL eneg_R = 0.5*dens_R*SQR(velx_R) + pres_R/(GAMMA-1.0);

    set_eigen_value(lambda, dens_L, velx_L, pres_L);
    set_LR_matrix(L, R, dens_L, velx_L, pres_L);
    calc_flux_jacobian_plus(Wp, L, R, lambda);

    set_eigen_value(lambda, dens_R, velx_R, pres_R);
    set_LR_matrix(L, R, dens_R, velx_R, pres_R);
    calc_flux_jacobian_minus(Wm, L, R, lambda);

#define WP(a,b) Wp[(a)+3*(b)]
#define WM(a,b) Wm[(a)+3*(b)]
    Ep[ix][0]
      = WP(0,0)*dens_L + WP(0,1)*momx_L + WP(0,2)*eneg_L;
    Ep[ix][1]
      = WP(1,0)*dens_L + WP(1,1)*momx_L + WP(1,2)*eneg_L;
    Ep[ix][2]
      = WP(2,0)*dens_L + WP(2,1)*momx_L + WP(2,2)*eneg_L;

    Em[ix][0]
      = WM(0,0)*dens_R + WM(0,1)*momx_R + WM(0,2)*eneg_R;
    Em[ix][1]
      = WM(1,0)*dens_R + WM(1,1)*momx_R + WM(1,2)*eneg_R;
    Em[ix][2]
      = WM(2,0)*dens_R + WM(2,1)*momx_R + WM(2,2)*eneg_R;
#undef WP
#undef WM

  }

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixp = ix+1;
    flux_iph[ix].dens_flux = Ep[ix][0] + Em[ixp][0];
    flux_iph[ix].momx_flux = Ep[ix][1] + Em[ixp][1];
    flux_iph[ix].eneg_flux = Ep[ix][2] + Em[ixp][2];
  }

  free(Ep);
  free(Em);
}

void calc_flux_FVS_MP5(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
		       struct run_param *this_run)
{
  static REAL L[9], R[9], Wp[9], Wm[9];
  static REAL wim2[3],wim1[3],wi[3],wip1[3],wip2[3];
  static REAL wiph[3],wimh[3];
  static REAL lambda[3];

  REAL (*Ep)[3], (*Em)[3];

  Ep = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  Em = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixm2 = ix-2; if(ixm2<0) ixm2 = 0;
    int ixm1 = ix-1; if(ixm1<0) ixm1 = 0;
    int ixp1 = ix+1; if(ixp1>this_run->nmesh-1) ixp1 = this_run->nmesh-1;
    int ixp2 = ix+2; if(ixp2>this_run->nmesh-1) ixp2 = this_run->nmesh-1;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);
    
    REAL dens_i = mesh[ix].dens;
    REAL velx_i = mesh[ix].momx/mesh[ix].dens;
    REAL eneg_i = mesh[ix].eneg;
    REAL pres_i = (GAMMA-1)*(mesh[ix].eneg - 0.5*dens_i*SQR(velx_i));

    set_LR_matrix(L, R, dens_i, velx_i, pres_i);

#define LL(a,b) L[(a)+3*(b)]  
    for(int m=0;m<3;m++) {
      wim2[m] = LL(m,0)*mesh[ixm2].dens 
	      + LL(m,1)*mesh[ixm2].momx 
	      + LL(m,2)*mesh[ixm2].eneg;
      wim1[m] = LL(m,0)*mesh[ixm1].dens 
	      + LL(m,1)*mesh[ixm1].momx 
	      + LL(m,2)*mesh[ixm1].eneg;
      wi[m]   = LL(m,0)*mesh[ix  ].dens 
	      + LL(m,1)*mesh[ix  ].momx 
	      + LL(m,2)*mesh[ix  ].eneg;
      wip1[m] = LL(m,0)*mesh[ixp1].dens 
	      + LL(m,1)*mesh[ixp1].momx 
	      + LL(m,2)*mesh[ixp1].eneg;
      wip2[m] = LL(m,0)*mesh[ixp2].dens 
	      + LL(m,1)*mesh[ixp2].momx 
	      + LL(m,2)*mesh[ixp2].eneg;
    }
#undef LL

    for(int m=0;m<3;m++) {
      wiph[m] = mp5_L(wim2[m], wim1[m], wi[m], wip1[m], wip2[m]);
      wimh[m] = mp5_R(wim2[m], wim1[m], wi[m], wip1[m], wip2[m]);
    }

#define RR(a,b) R[(a)+3*(b)]
    REAL dens_L = RR(0,0)*wiph[0] + RR(0,1)*wiph[1] + RR(0,2)*wiph[2];
    REAL momx_L = RR(1,0)*wiph[0] + RR(1,1)*wiph[1] + RR(1,2)*wiph[2];
    REAL eneg_L = RR(2,0)*wiph[0] + RR(2,1)*wiph[1] + RR(2,2)*wiph[2];
    REAL pres_L = (GAMMA-1.0)*(eneg_L - 0.5*SQR(momx_L)/dens_L);

    REAL dens_R = RR(0,0)*wimh[0] + RR(0,1)*wimh[1] + RR(0,2)*wimh[2];
    REAL momx_R = RR(1,0)*wimh[0] + RR(1,1)*wimh[1] + RR(1,2)*wimh[2];
    REAL eneg_R = RR(2,0)*wimh[0] + RR(2,1)*wimh[1] + RR(2,2)*wimh[2];
    REAL pres_R = (GAMMA-1.0)*(eneg_R - 0.5*SQR(momx_R)/dens_R);
#undef RR

    set_eigen_value(lambda, dens_L, momx_L/dens_L, pres_L);
    set_LR_matrix(L, R, dens_L, momx_L/dens_L, pres_L);
    calc_flux_jacobian_plus(Wp, L, R, lambda);

    set_eigen_value(lambda, dens_R, momx_R/dens_R, pres_R);
    set_LR_matrix(L, R, dens_R, momx_R/dens_R, pres_R);
    calc_flux_jacobian_minus(Wm, L, R, lambda);

#define WP(a,b) Wp[(a)+3*(b)]
#define WM(a,b) Wm[(a)+3*(b)]
    Ep[ix][0]
      = WP(0,0)*dens_L + WP(0,1)*momx_L + WP(0,2)*eneg_L;
    Ep[ix][1]
      = WP(1,0)*dens_L + WP(1,1)*momx_L + WP(1,2)*eneg_L;
    Ep[ix][2]
      = WP(2,0)*dens_L + WP(2,1)*momx_L + WP(2,2)*eneg_L;

    Em[ix][0]
      = WM(0,0)*dens_R + WM(0,1)*momx_R + WM(0,2)*eneg_R;
    Em[ix][1]
      = WM(1,0)*dens_R + WM(1,1)*momx_R + WM(1,2)*eneg_R;
    Em[ix][2]
      = WM(2,0)*dens_R + WM(2,1)*momx_R + WM(2,2)*eneg_R;
#undef WP
#undef WM

  }

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixp = ix+1;
    flux_iph[ix].dens_flux = Ep[ix][0] + Em[ixp][0];
    flux_iph[ix].momx_flux = Ep[ix][1] + Em[ixp][1];
    flux_iph[ix].eneg_flux = Ep[ix][2] + Em[ixp][2];
  }

  free(Ep);
  free(Em);
}

void calc_flux_FVS_MUSCL(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
			 struct run_param *this_run)
{
  static REAL L[9], R[9], Wp[9], Wm[9];
  static REAL lambda[3];

  REAL (*Ep)[3], (*Em)[3];

  Ep = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  Em = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixm1, ixp1;

    ixm1 = ix-1; if(ixm1 < 0) ixm1 = 0;
    ixp1 = ix+1; if(ixp1 > this_run->nmesh-1) ixp1 = this_run->nmesh - 1;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    REAL dens_im1, dens_i, dens_ip1;
    REAL velx_im1, velx_i, velx_ip1;
    REAL pres_im1, pres_i, pres_ip1;

    dens_im1 = mesh[ixm1].dens;
    dens_i   = mesh[ix  ].dens;
    dens_ip1 = mesh[ixp1].dens;

    velx_im1 = mesh[ixm1].momx/mesh[ixm1].dens;
    velx_i   = mesh[ix  ].momx/mesh[ix  ].dens;
    velx_ip1 = mesh[ixp1].momx/mesh[ixp1].dens;

    pres_im1 = (GAMMA-1.0)*(mesh[ixm1].eneg - 0.5*dens_im1*SQR(velx_im1));
    pres_i   = (GAMMA-1.0)*(mesh[ix  ].eneg - 0.5*dens_i  *SQR(velx_i  ));
    pres_ip1 = (GAMMA-1.0)*(mesh[ixp1].eneg - 0.5*dens_ip1*SQR(velx_ip1));

    REAL dens_L, dens_R;
    REAL velx_L, velx_R;
    REAL pres_L, pres_R;
    REAL eneg_L, eneg_R;

    dens_L = muscl_L(dens_im1, dens_i, dens_ip1);
    velx_L = muscl_L(velx_im1, velx_i, velx_ip1);
    pres_L = muscl_L(pres_im1, pres_i, pres_ip1);
    eneg_L = 0.5*SQR(velx_L)*dens_L + pres_L/(GAMMA-1.0);

    set_eigen_value(lambda, dens_L, velx_L, pres_L);
    set_LR_matrix(L, R, dens_L, velx_L, pres_L);
    calc_flux_jacobian_plus(Wp, L, R, lambda);

    dens_R = muscl_R(dens_im1, dens_i, dens_ip1);
    velx_R = muscl_R(velx_im1, velx_i, velx_ip1);
    pres_R = muscl_R(pres_im1, pres_i, pres_ip1);
    eneg_R = 0.5*SQR(velx_R)*dens_R + pres_R/(GAMMA-1.0);

    set_eigen_value(lambda, dens_R, velx_R, pres_R);
    set_LR_matrix(L, R, dens_R, velx_R, pres_R);
    calc_flux_jacobian_minus(Wm, L, R, lambda);

#define WP(a,b) Wp[(a)+3*(b)]
#define WM(a,b) Wm[(a)+3*(b)]
    Ep[ix][0]
      = WP(0,0)*dens_L + WP(0,1)*dens_L*velx_L + WP(0,2)*eneg_L;
    Ep[ix][1]
      = WP(1,0)*dens_L + WP(1,1)*dens_L*velx_L + WP(1,2)*eneg_L;
    Ep[ix][2]
      = WP(2,0)*dens_L + WP(2,1)*dens_L*velx_L + WP(2,2)*eneg_L;

    Em[ix][0]
      = WM(0,0)*dens_R + WM(0,1)*dens_R*velx_R + WM(0,2)*eneg_R;
    Em[ix][1]
      = WM(1,0)*dens_R + WM(1,1)*dens_R*velx_R + WM(1,2)*eneg_R;
    Em[ix][2]
      = WM(2,0)*dens_R + WM(2,1)*dens_R*velx_R + WM(2,2)*eneg_R;
#undef WP
#undef WM

  }

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixp = ix+1;
    flux_iph[ix].dens_flux = Ep[ix][0] + Em[ixp][0];
    flux_iph[ix].momx_flux = Ep[ix][1] + Em[ixp][1];
    flux_iph[ix].eneg_flux = Ep[ix][2] + Em[ixp][2];
  }

  free(Ep);
  free(Em);
}

void calc_flux_FVS_MUSCL2(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
			  struct run_param *this_run)
{
  static REAL L[9], R[9], Wp[9], Wm[9];
  static REAL lambda[3];

  REAL (*Ep)[3], (*Em)[3];

  Ep = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  Em = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixm1, ixp1;

    ixm1 = ix-1; if(ixm1 < 0) ixm1 = 0;
    ixp1 = ix+1; if(ixp1 > this_run->nmesh-1) ixp1 = this_run->nmesh - 1;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    REAL dens_im1, dens_i, dens_ip1;
    REAL momx_im1, momx_i, momx_ip1;
    REAL eneg_im1, eneg_i, eneg_ip1;

    dens_im1 = mesh[ixm1].dens;
    dens_i   = mesh[ix  ].dens;
    dens_ip1 = mesh[ixp1].dens;

    momx_im1 = mesh[ixm1].momx;
    momx_i   = mesh[ix  ].momx;
    momx_ip1 = mesh[ixp1].momx;

    eneg_im1 = mesh[ixm1].eneg;
    eneg_i   = mesh[ix  ].eneg;
    eneg_ip1 = mesh[ixp1].eneg;

    REAL dens_L, dens_R;
    REAL momx_L, momx_R;
    REAL velx_L, velx_R;
    REAL pres_L, pres_R;
    REAL eneg_L, eneg_R;

    dens_L = muscl_L(dens_im1, dens_i, dens_ip1);
    momx_L = muscl_L(momx_im1, momx_i, momx_ip1);
    velx_L = momx_L/dens_L;
    eneg_L = muscl_L(eneg_im1, eneg_i, eneg_ip1);
    pres_L = (eneg_L - 0.5*SQR(momx_L)/dens_L)*(GAMMA-1.0);

    set_eigen_value(lambda, dens_L, velx_L, pres_L);
    set_LR_matrix(L, R, dens_L, velx_L, pres_L);
    calc_flux_jacobian_plus(Wp, L, R, lambda);

    dens_R = muscl_R(dens_im1, dens_i, dens_ip1);
    momx_R = muscl_R(momx_im1, momx_i, momx_ip1);
    velx_R = momx_R/dens_R;
    eneg_R = muscl_R(eneg_im1, eneg_i, eneg_ip1);
    pres_R = (eneg_R - 0.5*SQR(momx_R)/dens_R)*(GAMMA-1.0);

    set_eigen_value(lambda, dens_R, velx_R, pres_R);
    set_LR_matrix(L, R, dens_R, velx_R, pres_R);
    calc_flux_jacobian_minus(Wm, L, R, lambda);

#define WP(a,b) Wp[(a)+3*(b)]
#define WM(a,b) Wm[(a)+3*(b)]
    Ep[ix][0]
      = WP(0,0)*dens_L + WP(0,1)*dens_L*velx_L + WP(0,2)*eneg_L;
    Ep[ix][1]
      = WP(1,0)*dens_L + WP(1,1)*dens_L*velx_L + WP(1,2)*eneg_L;
    Ep[ix][2]
      = WP(2,0)*dens_L + WP(2,1)*dens_L*velx_L + WP(2,2)*eneg_L;

    Em[ix][0]
      = WM(0,0)*dens_R + WM(0,1)*dens_R*velx_R + WM(0,2)*eneg_R;
    Em[ix][1]
      = WM(1,0)*dens_R + WM(1,1)*dens_R*velx_R + WM(1,2)*eneg_R;
    Em[ix][2]
      = WM(2,0)*dens_R + WM(2,1)*dens_R*velx_R + WM(2,2)*eneg_R;
#undef WP
#undef WM

  }

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixp = ix+1;
    flux_iph[ix].dens_flux = Ep[ix][0] + Em[ixp][0];
    flux_iph[ix].momx_flux = Ep[ix][1] + Em[ixp][1];
    flux_iph[ix].eneg_flux = Ep[ix][2] + Em[ixp][2];
  }

  free(Ep);
  free(Em);
}

void calc_flux_AUSMP_MP5(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
			 struct run_param *this_run)
{
  static REAL L[9], R[9], Wp[9], Wm[9];
  static REAL wim2[3],wim1[3],wi[3],wip1[3],wip2[3],wip3[3];
  static REAL wR[3],wL[3];
  static REAL lambda[3];

  REAL (*Ep)[3], (*Em)[3];

  Ep = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  Em = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixm2 = ix-2; if(ixm2<0) ixm2 = 0;
    int ixm1 = ix-1; if(ixm1<0) ixm1 = 0;
    int ixp1 = ix+1; if(ixp1>this_run->nmesh-1) ixp1 = this_run->nmesh-1;
    int ixp2 = ix+2; if(ixp2>this_run->nmesh-1) ixp2 = this_run->nmesh-1;
    int ixp3 = ix+3; if(ixp3>this_run->nmesh-1) ixp3 = this_run->nmesh-1;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);
    
    REAL dens_i = mesh[ix].dens;
    REAL velx_i = mesh[ix].momx/mesh[ix].dens;
    REAL eneg_i = mesh[ix].eneg;
    REAL pres_i = (GAMMA-1)*(mesh[ix].eneg - 0.5*dens_i*SQR(velx_i));

    set_LR_matrix(L, R, dens_i, velx_i, pres_i);

#define LL(a,b) L[(a)+3*(b)]  
    for(int m=0;m<3;m++) {
      wim2[m] = LL(m,0)*mesh[ixm2].dens 
	      + LL(m,1)*mesh[ixm2].momx 
	      + LL(m,2)*mesh[ixm2].eneg;
      wim1[m] = LL(m,0)*mesh[ixm1].dens 
	      + LL(m,1)*mesh[ixm1].momx 
	      + LL(m,2)*mesh[ixm1].eneg;
      wi[m]   = LL(m,0)*mesh[ix  ].dens 
	      + LL(m,1)*mesh[ix  ].momx 
	      + LL(m,2)*mesh[ix  ].eneg;
      wip1[m] = LL(m,0)*mesh[ixp1].dens 
	      + LL(m,1)*mesh[ixp1].momx 
	      + LL(m,2)*mesh[ixp1].eneg;
      wip2[m] = LL(m,0)*mesh[ixp2].dens 
	      + LL(m,1)*mesh[ixp2].momx 
	      + LL(m,2)*mesh[ixp2].eneg;
    }
#undef LL

    for(int m=0;m<3;m++) {
      wL[m] = mp5_L(wim2[m], wim1[m], wi[m], wip1[m], wip2[m]);
    }

#define RR(a,b) R[(a)+3*(b)]
    REAL dens_L = RR(0,0)*wL[0] + RR(0,1)*wL[1] + RR(0,2)*wL[2];
    REAL momx_L = RR(1,0)*wL[0] + RR(1,1)*wL[1] + RR(1,2)*wL[2];
    REAL eneg_L = RR(2,0)*wL[0] + RR(2,1)*wL[1] + RR(2,2)*wL[2];
    REAL velx_L = momx_L/dens_L;
    REAL pres_L = (GAMMA-1.0)*(eneg_L - 0.5*SQR(momx_L)/dens_L);
    REAL enth_L = (eneg_L + pres_L)/dens_L;
    REAL etrp_L = pres_L/pow(dens_L,GAMMA-1.0);
#undef RR

    REAL dens_ip1 = mesh[ixp1].dens;
    REAL velx_ip1 = mesh[ixp1].momx/mesh[ixp1].dens;
    REAL eneg_ip1 = mesh[ixp1].eneg;
    REAL pres_ip1 = (GAMMA-1)*(mesh[ixp1].eneg - 0.5*dens_ip1*SQR(velx_ip1));

    set_LR_matrix(L, R, dens_ip1, velx_ip1, pres_ip1);

#define LL(a,b) L[(a)+3*(b)]  
    for(int m=0;m<3;m++) {
      wim1[m] = LL(m,0)*mesh[ixm1].dens 
	      + LL(m,1)*mesh[ixm1].momx 
	      + LL(m,2)*mesh[ixm1].eneg;
      wi[m] = LL(m,0)*mesh[ix].dens 
	      + LL(m,1)*mesh[ix].momx 
	      + LL(m,2)*mesh[ix].eneg;
      wip1[m]   = LL(m,0)*mesh[ixp1].dens 
	      + LL(m,1)*mesh[ixp1].momx 
	      + LL(m,2)*mesh[ixp1].eneg;
      wip2[m] = LL(m,0)*mesh[ixp2].dens 
	      + LL(m,1)*mesh[ixp2].momx 
	      + LL(m,2)*mesh[ixp2].eneg;
      wip3[m] = LL(m,0)*mesh[ixp3].dens 
	      + LL(m,1)*mesh[ixp3].momx 
	      + LL(m,2)*mesh[ixp3].eneg;
    }
#undef LL

    for(int m=0;m<3;m++) {
      wR[m] = mp5_R(wim1[m], wi[m], wip1[m], wip2[m], wip3[m]);
    }

#define RR(a,b) R[(a)+3*(b)]
    REAL dens_R = RR(0,0)*wR[0] + RR(0,1)*wR[1] + RR(0,2)*wR[2];
    REAL momx_R = RR(1,0)*wR[0] + RR(1,1)*wR[1] + RR(1,2)*wR[2];
    REAL eneg_R = RR(2,0)*wR[0] + RR(2,1)*wR[1] + RR(2,2)*wR[2];
    REAL velx_R = momx_R/dens_R;
    REAL pres_R = (GAMMA-1.0)*(eneg_R - 0.5*SQR(momx_R)/dens_R);
    REAL enth_R = (eneg_R + pres_R)/dens_R;
    REAL etrp_R = pres_R/pow(dens_R,GAMMA-1.0);
#undef RR

    REAL cs;

    cs = sqrt(2.0*(GAMMA-1.0)/(GAMMA+1.0)*enth_L);
    REAL cs_L = SQR(cs)/fmax(cs, fabs(velx_L));    
    cs = sqrt(2.0*(GAMMA-1.0)/(GAMMA+1.0)*enth_R);
    REAL cs_R = SQR(cs)/fmax(cs, fabs(velx_R));

    REAL cs_h = fmin(cs_L, cs_R);
    REAL Mach_L = velx_L/cs_h;
    REAL Mach_R = velx_R/cs_h;

    REAL Mach_iph = M_p(Mach_L) + M_m(Mach_R);
    REAL pres_iph = P_p(Mach_L)*pres_L + P_m(Mach_R)*pres_R;

    if(Mach_iph > 0.0) {
      flux_iph[ix].dens_flux = cs_h*Mach_iph*dens_L;
      flux_iph[ix].momx_flux = cs_h*Mach_iph*dens_L*velx_L + pres_iph;
      flux_iph[ix].eneg_flux = cs_h*Mach_iph*dens_L*enth_L;
      flux_iph[ix].etrp_flux = cs_h*Mach_iph*etrp_L;
    }else{
      flux_iph[ix].dens_flux = cs_h*Mach_iph*dens_R;
      flux_iph[ix].momx_flux = cs_h*Mach_iph*dens_R*velx_R + pres_iph;
      flux_iph[ix].eneg_flux = cs_h*Mach_iph*dens_R*enth_R;
      flux_iph[ix].etrp_flux = cs_h*Mach_iph*etrp_R;
    }
  }
}

void calc_flux_AUSMP_MUSCL(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
			   struct run_param *this_run)
{
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL cs_h, cs;
    int ixm1, ixp1, ixp2;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    ixm1 = ix-1; if(ixm1 < 0) ixm1 = 0;
    ixp1 = ix+1; if(ixp1 > this_run->nmesh-1) ixp1 = this_run->nmesh - 1;
    ixp2 = ix+2; if(ixp2 > this_run->nmesh-1) ixp2 = this_run->nmesh - 1;

    REAL dens_im1, dens_i, dens_ip1, dens_ip2;
    REAL velx_im1, velx_i, velx_ip1, velx_ip2;
    REAL pres_im1, pres_i, pres_ip1, pres_ip2;
    REAL enth_i, enth_ip1;

    dens_im1 = mesh[ixm1].dens;
    dens_i   = mesh[ix  ].dens;
    dens_ip1 = mesh[ixp1].dens;
    dens_ip2 = mesh[ixp2].dens;

    velx_im1 = mesh[ixm1].momx/mesh[ixm1].dens;
    velx_i   = mesh[ix  ].momx/mesh[ix  ].dens;
    velx_ip1 = mesh[ixp1].momx/mesh[ixp1].dens;
    velx_ip2 = mesh[ixp2].momx/mesh[ixp2].dens;

    pres_im1 = (GAMMA-1.0)*(mesh[ixm1].eneg - 0.5*dens_im1*SQR(velx_im1));
    pres_i   = (GAMMA-1.0)*(mesh[ix  ].eneg - 0.5*dens_i  *SQR(velx_i  ));
    pres_ip1 = (GAMMA-1.0)*(mesh[ixp1].eneg - 0.5*dens_ip1*SQR(velx_ip1));
    pres_ip2 = (GAMMA-1.0)*(mesh[ixp2].eneg - 0.5*dens_ip2*SQR(velx_ip2));

    enth_i   = (mesh[ix].eneg + pres_i)/dens_i;
    enth_ip1 = (mesh[ixp1].eneg + pres_ip1)/dens_ip1;

    REAL dens_L, dens_R;
    REAL velx_L, velx_R;
    REAL pres_L, pres_R;
    REAL eneg_L, eneg_R;
    REAL enth_L, enth_R;
    REAL etrp_L, etrp_R;

    dens_L = muscl_L(dens_im1, dens_i, dens_ip1);
    velx_L = muscl_L(velx_im1, velx_i, velx_ip1);
    pres_L = muscl_L(pres_im1, pres_i, pres_ip1);
    eneg_L = 0.5*dens_L*SQR(velx_L) + pres_L/(GAMMA-1.0);
    enth_L = (eneg_L + pres_L)/dens_L;
    etrp_L = pres_L/pow(dens_L,GAMMA-1.0);

    dens_R = muscl_R(dens_i, dens_ip1, dens_ip2);
    velx_R = muscl_R(velx_i, velx_ip1, velx_ip2);
    pres_R = muscl_R(pres_i, pres_ip1, pres_ip2);
    eneg_R = 0.5*dens_R*SQR(velx_R) + pres_R/(GAMMA-1.0);
    enth_R = (eneg_R + pres_R)/dens_R;
    etrp_R = pres_R/pow(dens_R,GAMMA-1.0);

    cs = sqrt(2.0*(GAMMA-1.0)/(GAMMA+1.0)*enth_L);
    REAL cs_L = SQR(cs)/fmax(cs, fabs(velx_L));    
    cs = sqrt(2.0*(GAMMA-1.0)/(GAMMA+1.0)*enth_R);
    REAL cs_R = SQR(cs)/fmax(cs, fabs(velx_R));

    cs_h = fmin(cs_L, cs_R);
    REAL Mach_L = velx_L/cs_h;
    REAL Mach_R = velx_R/cs_h;

    REAL Mach_iph = M_p(Mach_L) + M_m(Mach_R);
    REAL pres_iph = P_p(Mach_L)*pres_L + P_m(Mach_R)*pres_R;

    if(Mach_iph > 0.0) {
      flux_iph[ix].dens_flux = cs_h*Mach_iph*dens_L;
      flux_iph[ix].momx_flux = cs_h*Mach_iph*dens_L*velx_L + pres_iph;
      flux_iph[ix].eneg_flux = cs_h*Mach_iph*dens_L*enth_L;
      flux_iph[ix].etrp_flux = cs_h*Mach_iph*etrp_L;
    }else{
      flux_iph[ix].dens_flux = cs_h*Mach_iph*dens_R;
      flux_iph[ix].momx_flux = cs_h*Mach_iph*dens_R*velx_R + pres_iph;
      flux_iph[ix].eneg_flux = cs_h*Mach_iph*dens_R*enth_R;
      flux_iph[ix].etrp_flux = cs_h*Mach_iph*etrp_R;
    }

  }
}

void calc_flux_AUSMP(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
		     struct run_param *this_run)
{
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL cs_h,cs;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    /* compute the flux at i+1/2*/
    REAL dens_L = mesh[ix].dens;
    REAL velx_L = mesh[ix].momx/dens_L;
    REAL pres_L = (GAMMA-1.0)*(mesh[ix].eneg - 0.5*dens_L*SQR(velx_L));
    REAL enth_L = (mesh[ix].eneg + pres_L)/dens_L;
    REAL etrp_L = pres_L/pow(dens_L,GAMMA-1.0);
    
    cs = sqrt(2.0*(GAMMA-1.0)/(GAMMA+1.0)*enth_L);
    REAL cs_L = SQR(cs)/fmax(cs, fabs(velx_L));

    int ixp = ix + 1; if(ixp>this_run->nmesh-1) ixp=this_run->nmesh-1;
    REAL dens_R = mesh[ixp].dens;
    REAL velx_R = mesh[ixp].momx/dens_R;
    REAL pres_R = (GAMMA-1.0)*(mesh[ixp].eneg - 0.5*dens_R*SQR(velx_R));
    REAL enth_R = (mesh[ixp].eneg + pres_R)/dens_R;
    REAL etrp_R = pres_R/pow(dens_R,GAMMA-1.0);

    cs = sqrt(2.0*(GAMMA-1.0)/(GAMMA+1.0)*enth_R);
    REAL cs_R = SQR(cs)/fmax(cs, fabs(velx_R));

    cs_h = fmin(cs_L, cs_R);
    REAL Mach_L = velx_L/cs_h;
    REAL Mach_R = velx_R/cs_h;

    REAL Mach_iph = M_p(Mach_L) + M_m(Mach_R);
    REAL pres_iph = P_p(Mach_L)*pres_L + P_m(Mach_R)*pres_R;

    if(Mach_iph >= 0.0) {
      flux_iph[ix].dens_flux = cs_h*Mach_iph*dens_L;
      flux_iph[ix].momx_flux = cs_h*Mach_iph*dens_L*velx_L + pres_iph;
      flux_iph[ix].eneg_flux = cs_h*Mach_iph*dens_L*enth_L;
      flux_iph[ix].etrp_flux = cs_h*Mach_iph*etrp_L;
    }else{
      flux_iph[ix].dens_flux = cs_h*Mach_iph*dens_R;
      flux_iph[ix].momx_flux = cs_h*Mach_iph*dens_R*velx_R + pres_iph;
      flux_iph[ix].eneg_flux = cs_h*Mach_iph*dens_R*enth_R;
      flux_iph[ix].etrp_flux = cs_h*Mach_iph*dens_R;
    }
  }
}

void calc_flux_FVS(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
		   struct run_param *this_run)
{
  static REAL L[9], R[9], Wp[9], Wm[9];
  static REAL lambda[3];

  REAL (*Ep)[3], (*Em)[3];

  Ep = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  Em = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);

  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL velx = mesh[ix].momx/mesh[ix].dens;
    REAL pres = (GAMMA-1.0)*(mesh[ix].eneg-0.5*SQR(velx)*mesh[ix].dens);

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    set_eigen_value(lambda, mesh[ix].dens, velx, pres);

    set_LR_matrix(L, R, mesh[ix].dens, velx, pres);
    calc_flux_jacobian(Wp, Wm, L, R, lambda);

#define WP(a,b) Wp[(a)+3*(b)]
#define WM(a,b) Wm[(a)+3*(b)]
    Ep[ix][0]
      = WP(0,0)*mesh[ix].dens + WP(0,1)*mesh[ix].momx + WP(0,2)*mesh[ix].eneg;
    Ep[ix][1]
      = WP(1,0)*mesh[ix].dens + WP(1,1)*mesh[ix].momx + WP(1,2)*mesh[ix].eneg;
    Ep[ix][2]
      = WP(2,0)*mesh[ix].dens + WP(2,1)*mesh[ix].momx + WP(2,2)*mesh[ix].eneg;

    Em[ix][0]
      = WM(0,0)*mesh[ix].dens + WM(0,1)*mesh[ix].momx + WM(0,2)*mesh[ix].eneg;
    Em[ix][1]
      = WM(1,0)*mesh[ix].dens + WM(1,1)*mesh[ix].momx + WM(1,2)*mesh[ix].eneg;
    Em[ix][2]
      = WM(2,0)*mesh[ix].dens + WM(2,1)*mesh[ix].momx + WM(2,2)*mesh[ix].eneg;
#undef WP
#undef WM
  }

#ifdef __PERIODIC__
  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixp = ix+1;if(ixp>this_run->nmesh-1) ixp = 0;
    flux_iph[ix].dens_flux = Ep[ix][0] + Em[ixp][0];
    flux_iph[ix].momx_flux = Ep[ix][1] + Em[ixp][1];
    flux_iph[ix].eneg_flux = Ep[ix][2] + Em[ixp][2];
  }
#else
  for(int ix=0;ix<this_run->nmesh-1;ix++) {
    int ixp = ix+1;
    flux_iph[ix].dens_flux = Ep[ix][0] + Em[ixp][0];
    flux_iph[ix].momx_flux = Ep[ix][1] + Em[ixp][1];
    flux_iph[ix].eneg_flux = Ep[ix][2] + Em[ixp][2];
  }
#endif

  free(Ep);
  free(Em);
}

void calc_flux_HLLC_MP5(struct fluid_1d *mesh, struct fluid_flux *flux, 
			struct run_param *this_run)
{
  static REAL L[9], R[9], Wp[9], Wm[9];
  static REAL wim2[3],wim1[3],wi[3],wip1[3],wip2[3],wip3[3];
  static REAL wR[3],wL[3];
  static REAL lambda[3];

  REAL (*Ep)[3], (*Em)[3];

  Ep = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  Em = (REAL (*)[3]) malloc(sizeof(REAL)*this_run->nmesh*3);
  
  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixp3 = ix+3; if(ixp3>this_run->nmesh-1) ixp3 = this_run->nmesh-1;
    int ixp2 = ix+2; if(ixp2>this_run->nmesh-1) ixp2 = this_run->nmesh-1;
    int ixp1 = ix+1; if(ixp1>this_run->nmesh-1) ixp1 = this_run->nmesh-1;
    int ixm1 = ix-1; if(ixm1<0) ixm1 = 0;
    int ixm2 = ix-2; if(ixm2<0) ixm2 = 0;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    REAL dens_i = mesh[ix].dens;
    REAL velx_i = mesh[ix].momx/mesh[ix].dens;
    REAL eneg_i = mesh[ix].eneg;
    REAL pres_i = (GAMMA-1)*(mesh[ix].eneg - 0.5*dens_i*SQR(velx_i));

    set_LR_matrix(L, R, dens_i, velx_i, pres_i);

#define LL(a,b) L[(a)+3*(b)]  
    for(int m=0;m<3;m++) {
      wim2[m] = LL(m,0)*mesh[ixm2].dens 
	      + LL(m,1)*mesh[ixm2].momx 
	      + LL(m,2)*mesh[ixm2].eneg;
      wim1[m] = LL(m,0)*mesh[ixm1].dens 
	      + LL(m,1)*mesh[ixm1].momx 
	      + LL(m,2)*mesh[ixm1].eneg;
      wi[m]   = LL(m,0)*mesh[ix  ].dens 
	      + LL(m,1)*mesh[ix  ].momx 
	      + LL(m,2)*mesh[ix  ].eneg;
      wip1[m] = LL(m,0)*mesh[ixp1].dens 
	      + LL(m,1)*mesh[ixp1].momx 
	      + LL(m,2)*mesh[ixp1].eneg;
      wip2[m] = LL(m,0)*mesh[ixp2].dens 
	      + LL(m,1)*mesh[ixp2].momx 
	      + LL(m,2)*mesh[ixp2].eneg;
    }
#undef LL

    for(int m=0;m<3;m++) {
      wL[m] = mp5_L(wim2[m], wim1[m], wi[m], wip1[m], wip2[m]);
    }

#define RR(a,b) R[(a)+3*(b)]
    REAL dens_L = RR(0,0)*wL[0] + RR(0,1)*wL[1] + RR(0,2)*wL[2];
    REAL momx_L = RR(1,0)*wL[0] + RR(1,1)*wL[1] + RR(1,2)*wL[2];
    REAL eneg_L = RR(2,0)*wL[0] + RR(2,1)*wL[1] + RR(2,2)*wL[2];
    REAL velx_L = momx_L/dens_L;
    REAL pres_L = (GAMMA-1.0)*(eneg_L - 0.5*SQR(momx_L)/dens_L);
    REAL enth_L = (eneg_L + pres_L)/dens_L;
    REAL etrp_L = pres_L/pow(dens_L,GAMMA-1.0);
#undef RR

    REAL dens_ip1 = mesh[ixp1].dens;
    REAL velx_ip1 = mesh[ixp1].momx/mesh[ixp1].dens;
    REAL eneg_ip1 = mesh[ixp1].eneg;
    REAL pres_ip1 = (GAMMA-1)*(mesh[ixp1].eneg - 0.5*dens_ip1*SQR(velx_ip1));

    set_LR_matrix(L, R, dens_ip1, velx_ip1, pres_ip1);

#define LL(a,b) L[(a)+3*(b)]  
    for(int m=0;m<3;m++) {
      wim1[m] = LL(m,0)*mesh[ixm1].dens 
	      + LL(m,1)*mesh[ixm1].momx 
	      + LL(m,2)*mesh[ixm1].eneg;
      wi[m] = LL(m,0)*mesh[ix].dens 
	      + LL(m,1)*mesh[ix].momx 
	      + LL(m,2)*mesh[ix].eneg;
      wip1[m]   = LL(m,0)*mesh[ixp1].dens 
	      + LL(m,1)*mesh[ixp1].momx 
	      + LL(m,2)*mesh[ixp1].eneg;
      wip2[m] = LL(m,0)*mesh[ixp2].dens 
	      + LL(m,1)*mesh[ixp2].momx 
	      + LL(m,2)*mesh[ixp2].eneg;
      wip3[m] = LL(m,0)*mesh[ixp3].dens 
	      + LL(m,1)*mesh[ixp3].momx 
	      + LL(m,2)*mesh[ixp3].eneg;
    }
#undef LL

    for(int m=0;m<3;m++) {
      wR[m] = mp5_R(wim1[m], wi[m], wip1[m], wip2[m], wip3[m]);
    }

#define RR(a,b) R[(a)+3*(b)]
    REAL dens_R = RR(0,0)*wR[0] + RR(0,1)*wR[1] + RR(0,2)*wR[2];
    REAL momx_R = RR(1,0)*wR[0] + RR(1,1)*wR[1] + RR(1,2)*wR[2];
    REAL eneg_R = RR(2,0)*wR[0] + RR(2,1)*wR[1] + RR(2,2)*wR[2];
    REAL velx_R = momx_R/dens_R;
    REAL pres_R = (GAMMA-1.0)*(eneg_R - 0.5*SQR(momx_R)/dens_R);
    REAL enth_R = (eneg_R + pres_R)/dens_R;
    REAL etrp_R = pres_R/pow(dens_R,GAMMA-1.0);
#undef RR

    REAL cs_L = sqrt(GAMMA*pres_L/dens_L);
    REAL cs_R = sqrt(GAMMA*pres_L/dens_R);
    
    /* pressure estimate */
    REAL pres_star;
#if 1
    /* acoustic type arpproximation */
    REAL pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_star = fmax(0.0, pres_pvrs);
#endif
#if 1
    /* two-rarefaction riemann solver */
    REAL zindx = 0.5*(GAMMA-1.0)/GAMMA;
    pres_star = pow((cs_L+cs_R-0.5*(GAMMA-1.0)*(velx_R-velx_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#if 0
    /* two-shock riemann solver */
    REAL pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    REAL pres_0 = fmax(0.0, pres_pvrs);
    REAL gL, gR;
    gL = sqrt(2.0/((GAMMA+1.0)*dens_L)/(pres_0+(GAMMA-1.0)/(GAMMA+1.0)*pres_L));
    gR = sqrt(2.0/((GAMMA+1.0)*dens_R)/(pres_0+(GAMMA-1.0)/(GAMMA+1.0)*pres_R));

    pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif

    REAL qL, qR;

    if(pres_star <= pres_L) {
      qL = 1.0;
    }else{
      qL = sqrt(1.0+0.5*(GAMMA+1.0)/GAMMA*(pres_star/pres_L-1.0));
    }
    
    if(pres_star <= pres_R) {
      qR = 1.0;
    }else{
      qR = sqrt(1.0+0.5*(GAMMA+1.0)/GAMMA*(pres_star/pres_R-1.0));
    }

    REAL SL, Sstar, SR;
    SL = fmin(velx_L - cs_L*qL, velx_R - cs_R*qR);
    SR = fmax(velx_R + cs_R*qR, velx_L + cs_L*qL);

    Sstar = (pres_R-pres_L + dens_L*velx_L*(SL-velx_L) - dens_R*velx_R*(SR-velx_R))/(dens_L*(SL-velx_L)-dens_R*(SR-velx_R));

    if(0.0<=SL){
      flux[ix].dens_flux = dens_L*velx_L;
      flux[ix].momx_flux = pres_L + dens_L*SQR(velx_L);
      flux[ix].eneg_flux = (eneg_L+pres_L)*velx_L;
      flux[ix].etrp_flux = etrp_L*velx_L;
    }else if(SL<0.0 && 0.0<=Sstar) {
      REAL pres_L_star = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
      REAL pres_R_star = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
      REAL pres_LR_star = 0.5*(pres_L_star + pres_R_star);
      flux[ix].dens_flux = Sstar*(SL*dens_L - dens_L*velx_L)/(SL-Sstar);
      flux[ix].momx_flux = 
	(Sstar*(SL*momx_L - (pres_L+dens_L*SQR(velx_L))) + SL*pres_L_star)/(SL-Sstar);
      flux[ix].eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*velx_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
      flux[ix].etrp_flux = Sstar*(SL*etrp_L - etrp_L*velx_L)/(SL-Sstar);
    }else if(Sstar < 0.0 && 0.0 <=SR) {
      REAL pres_L_star = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
      REAL pres_R_star = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
      REAL pres_LR_star = 0.5*(pres_L_star + pres_R_star);
      flux[ix].dens_flux = Sstar*(SR*dens_R - dens_R*velx_R)/(SR-Sstar);
      flux[ix].momx_flux = (Sstar*(SR*momx_R - (pres_R+dens_R*SQR(velx_R))) + SR*pres_R_star)/(SR-Sstar);
      flux[ix].eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*velx_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
      flux[ix].etrp_flux = Sstar*(SR*etrp_R - etrp_R*velx_R)/(SR-Sstar);
    }else if(SR<0.0) {
      flux[ix].dens_flux = dens_R*velx_R;
      flux[ix].momx_flux = pres_R + dens_R*SQR(velx_R);
      flux[ix].eneg_flux = (eneg_R+pres_R)*velx_R;
      flux[ix].etrp_flux = etrp_R*velx_R;
    }

  }

}

void calc_flux_HLLC_MUSCL(struct fluid_1d *mesh, struct fluid_flux *flux, 
			  struct run_param *this_run)
{
  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixp2 = ix+2; if(ixp2>this_run->nmesh-1) ixp2 = this_run->nmesh-1;
    int ixp1 = ix+1; if(ixp1>this_run->nmesh-1) ixp1 = this_run->nmesh-1;
    int ixm1 = ix-1; if(ixm1<0) ixm1 = 0;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    REAL dens_im1 = mesh[ixm1].dens;
    REAL dens_i   = mesh[ix].dens;
    REAL dens_ip1 = mesh[ixp1].dens;
    REAL dens_ip2 = mesh[ixp2].dens;
    REAL velx_im1 = mesh[ixm1].momx/mesh[ixm1].dens;
    REAL velx_i   = mesh[ix].momx/mesh[ix].dens;
    REAL velx_ip1 = mesh[ixp1].momx/mesh[ixp1].dens;
    REAL velx_ip2 = mesh[ixp2].momx/mesh[ixp2].dens;
    REAL pres_im1 = (GAMMA-1.0)*(mesh[ixm1].eneg-0.5*dens_im1*SQR(velx_im1));
    REAL pres_i   = (GAMMA-1.0)*(mesh[ix  ].eneg-0.5*dens_i  *SQR(velx_i  ));
    REAL pres_ip1 = (GAMMA-1.0)*(mesh[ixp1].eneg-0.5*dens_ip1*SQR(velx_ip1));
    REAL pres_ip2 = (GAMMA-1.0)*(mesh[ixp2].eneg-0.5*dens_ip2*SQR(velx_ip2));
    
    REAL dens_L = muscl_L(dens_im1, dens_i, dens_ip1);
    REAL velx_L = muscl_L(velx_im1, velx_i, velx_ip1);
    REAL pres_L = muscl_L(pres_im1, pres_i, pres_ip1);
    REAL eneg_L = 0.5*dens_L*SQR(velx_L) + pres_L/(GAMMA-1.0);
    REAL enth_L = (eneg_L + pres_L)/dens_L;
    REAL momx_L = dens_L*velx_L;
    REAL etrp_L = pres_L/pow(dens_L,GAMMA-1.0);

    REAL dens_R = muscl_R(dens_i, dens_ip1, dens_ip2);
    REAL velx_R = muscl_R(velx_i, velx_ip1, velx_ip2);
    REAL pres_R = muscl_R(pres_i, pres_ip1, pres_ip2);
    REAL eneg_R = 0.5*dens_R*SQR(velx_R) + pres_R/(GAMMA-1.0);
    REAL enth_R = (eneg_R + pres_R)/dens_R;
    REAL momx_R = dens_R*velx_R;
    REAL etrp_R = pres_R/pow(dens_R,GAMMA-1.0);

    REAL cs_L = sqrt(GAMMA*pres_L/dens_L);
    REAL cs_R = sqrt(GAMMA*pres_R/dens_R);
    
    /* pressure estimate */
    REAL pres_star;
#if 0
    /* acoustic type arpproximation */
    REAL pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_star = fmax(0.0, pres_pvrs);
#endif
#if 1
    /* two-rarefaction riemann solver */
    REAL zindx = 0.5*(GAMMA-1.0)/GAMMA;
    pres_star = pow((cs_L+cs_R-0.5*(GAMMA-1.0)*(velx_R-velx_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#if 0
    /* two-shock riemann solver */
    REAL pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    REAL pres_0 = fmax(0.0, pres_pvrs);
    REAL gL, gR;
    gL = sqrt(2.0/((GAMMA+1.0)*dens_L)/(pres_0+(GAMMA-1.0)/(GAMMA+1.0)*pres_L));
    gR = sqrt(2.0/((GAMMA+1.0)*dens_R)/(pres_0+(GAMMA-1.0)/(GAMMA+1.0)*pres_R));

    pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif

    REAL qL, qR;

    if(pres_star <= pres_L) {
      qL = 1.0;
    }else{
      qL = sqrt(1.0+0.5*(GAMMA+1.0)/GAMMA*(pres_star/pres_L-1.0));
    }
    
    if(pres_star <= pres_R) {
      qR = 1.0;
    }else{
      qR = sqrt(1.0+0.5*(GAMMA+1.0)/GAMMA*(pres_star/pres_R-1.0));
    }

    REAL SL, Sstar, SR;
    SL = fmin(velx_L - cs_L*qL, velx_R - cs_R*qR);
    SR = fmax(velx_R + cs_R*qR, velx_L + cs_L*qL);

    Sstar = (pres_R-pres_L + dens_L*velx_L*(SL-velx_L) - dens_R*velx_R*(SR-velx_R))/(dens_L*(SL-velx_L)-dens_R*(SR-velx_R));

    if(0.0<=SL){
      flux[ix].dens_flux = dens_L*velx_L;
      flux[ix].momx_flux = pres_L + dens_L*SQR(velx_L);
      flux[ix].eneg_flux = (eneg_L+pres_L)*velx_L;
      flux[ix].etrp_flux = etrp_L*velx_L;
    }else if(SL<0.0 && 0.0<=Sstar) {
      REAL pres_L_star = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
      REAL pres_R_star = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
      REAL pres_LR_star = 0.5*(pres_L_star + pres_R_star);
      flux[ix].dens_flux = Sstar*(SL*dens_L - dens_L*velx_L)/(SL-Sstar);
      flux[ix].momx_flux = 
	(Sstar*(SL*momx_L - (pres_L+dens_L*SQR(velx_L))) + SL*pres_L_star)/(SL-Sstar);
      flux[ix].eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*velx_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
      flux[ix].etrp_flux = Sstar*(SL*etrp_L - etrp_L*velx_L)/(SL-Sstar);
    }else if(Sstar < 0.0 && 0.0 <=SR) {
      REAL pres_L_star = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
      REAL pres_R_star = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
      REAL pres_LR_star = 0.5*(pres_L_star + pres_R_star);
      flux[ix].dens_flux = Sstar*(SR*dens_R - dens_R*velx_R)/(SR-Sstar);
      flux[ix].momx_flux = (Sstar*(SR*momx_R - (pres_R+dens_R*SQR(velx_R))) + SR*pres_R_star)/(SR-Sstar);
      flux[ix].eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*velx_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
      flux[ix].etrp_flux = Sstar*(SR*etrp_R - etrp_R*velx_R)/(SR-Sstar);
    }else if(SR<0.0) {
      flux[ix].dens_flux = dens_R*velx_R;
      flux[ix].momx_flux = pres_R + dens_R*SQR(velx_R);
      flux[ix].eneg_flux = (eneg_R+pres_R)*velx_R;
      flux[ix].etrp_flux = etrp_R*velx_R;
    }
  }

}

void calc_flux_HLLC(struct fluid_1d *mesh, struct fluid_flux *flux, 
		    struct run_param *this_run)
{
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL pres_L, pres_R;
    REAL dens_L, dens_R;
    REAL velx_L, velx_R;
    REAL eneg_L, eneg_R;
    REAL etrp_L, etrp_R;
    REAL momx_L, momx_R;
    REAL cs_L, cs_R;

    int ixp=ix+1; if(ixp>this_run->nmesh-1) ixp=this_run->nmesh-1;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    momx_L = mesh[ix].momx;
    velx_L = mesh[ix].momx/mesh[ix].dens;
    dens_L = mesh[ix].dens;
    pres_L = (mesh[ix].eneg-0.5*SQR(mesh[ix].momx)/mesh[ix].dens)*(GAMMA-1.0);
    eneg_L = mesh[ix].eneg;
    etrp_L = pres_L/pow(dens_L,GAMMA-1.0);

    momx_R = mesh[ixp].momx;
    velx_R = mesh[ixp].momx/mesh[ixp].dens;
    dens_R = mesh[ixp].dens;
    pres_R = (mesh[ixp].eneg-0.5*SQR(mesh[ixp].momx)/mesh[ixp].dens)*(GAMMA-1.0);
    eneg_R = mesh[ixp].eneg;
    etrp_R = pres_R/pow(dens_R,GAMMA-1.0);

    cs_L = sqrt(GAMMA*pres_L/dens_L);
    cs_R = sqrt(GAMMA*pres_R/dens_R);
    
    /* pressure estimate */
    REAL pres_star;
#if 0
    /* acoustic type arpproximation */
    REAL pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_star = fmax(0.0, pres_pvrs);
#endif
#if 1
    /* two-rarefaction riemann solver */
    REAL zindx = 0.5*(GAMMA-1.0)/GAMMA;
    pres_star = pow((cs_L+cs_R-0.5*(GAMMA-1.0)*(velx_R-velx_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#if 0
    /* two-shock riemann solver */
    REAL pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    REAL pres_0 = fmax(0.0, pres_pvrs);
    REAL gL, gR;
    gL = sqrt(2.0/((GAMMA+1.0)*dens_L)/(pres_0+(GAMMA-1.0)/(GAMMA+1.0)*pres_L));
    gR = sqrt(2.0/((GAMMA+1.0)*dens_R)/(pres_0+(GAMMA-1.0)/(GAMMA+1.0)*pres_R));

    pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif

    REAL qL, qR;

    if(pres_star <= pres_L) {
      qL = 1.0;
    }else{
      qL = sqrt(1.0+0.5*(GAMMA+1.0)/GAMMA*(pres_star/pres_L-1.0));
    }
    
    if(pres_star <= pres_R) {
      qR = 1.0;
    }else{
      qR = sqrt(1.0+0.5*(GAMMA+1.0)/GAMMA*(pres_star/pres_R-1.0));
    }

    REAL SL, Sstar, SR;
    SL = fmin(velx_L - cs_L*qL, velx_R - cs_R*qR);
    SR = fmax(velx_R + cs_R*qR, velx_L + cs_L*qL);

    Sstar = (pres_R-pres_L + dens_L*velx_L*(SL-velx_L) - dens_R*velx_R*(SR-velx_R))/(dens_L*(SL-velx_L)-dens_R*(SR-velx_R));

    if(0.0<=SL){
      flux[ix].dens_flux = dens_L*velx_L;
      flux[ix].momx_flux = pres_L + dens_L*SQR(velx_L);
      flux[ix].eneg_flux = (eneg_L+pres_L)*velx_L;
      flux[ix].etrp_flux = etrp_L*velx_L;
    }else if(SL<0.0 && 0.0<=Sstar) {
      REAL pres_L_star = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
      REAL pres_R_star = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
      REAL pres_LR_star = 0.5*(pres_L_star + pres_R_star);
      flux[ix].dens_flux = Sstar*(SL*dens_L - dens_L*velx_L)/(SL-Sstar);
      flux[ix].momx_flux = 
	(Sstar*(SL*momx_L - (pres_L+dens_L*SQR(velx_L))) + SL*pres_L_star)/(SL-Sstar);
      flux[ix].eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*velx_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
      flux[ix].etrp_flux = Sstar*(SL*etrp_L - etrp_L*velx_L)/(SL-Sstar);
    }else if(Sstar < 0.0 && 0.0 <=SR) {
      REAL pres_L_star = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
      REAL pres_R_star = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
      REAL pres_LR_star = 0.5*(pres_L_star + pres_R_star);
      flux[ix].dens_flux = Sstar*(SR*dens_R - dens_R*velx_R)/(SR-Sstar);
      flux[ix].momx_flux = (Sstar*(SR*momx_R - (pres_R+dens_R*SQR(velx_R))) + SR*pres_R_star)/(SR-Sstar);
      flux[ix].eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*velx_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
      flux[ix].etrp_flux = Sstar*(SR*etrp_R - etrp_R*velx_R)/(SR-Sstar);
    }else if(SR<0.0) {
      flux[ix].dens_flux = dens_R*velx_R;
      flux[ix].momx_flux = pres_R + dens_R*SQR(velx_R);
      flux[ix].eneg_flux = (eneg_R+pres_R)*velx_R;
      flux[ix].etrp_flux = etrp_R*velx_R;
    }
  }

}

void calc_flux_HLL(struct fluid_1d *mesh, struct fluid_flux *flux, 
		   struct run_param *this_run)
{
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL pres_L, pres_R;
    REAL dens_L, dens_R;
    REAL velx_L, velx_R;
    REAL eneg_L, eneg_R;
    REAL enth_L, enth_R;
    REAL etrp_L, etrp_R;
    REAL momx_L, momx_R;
    REAL cs_L, cs_R;

    int ixp=ix+1; if(ixp>this_run->nmesh-1) ixp=this_run->nmesh-1;

    mesh[ix].high_mach = high_mach_detector(&mesh[ix]);
    mesh[ix].under_shock = shock_detector(mesh, ix, this_run);

    momx_L = mesh[ix].momx;
    velx_L = mesh[ix].momx/mesh[ix].dens;
    dens_L = mesh[ix].dens;
    pres_L = (mesh[ix].eneg-0.5*SQR(mesh[ix].momx)/mesh[ix].dens)*(GAMMA-1.0);
    eneg_L = mesh[ix].eneg;
    enth_L = (eneg_L + pres_L)/dens_L;
    etrp_L = pres_L/pow(dens_L,GAMMA-1.0);

    momx_R = mesh[ixp].momx;
    velx_R = mesh[ixp].momx/mesh[ixp].dens;
    dens_R = mesh[ixp].dens;
    pres_R = (mesh[ixp].eneg-0.5*SQR(mesh[ixp].momx)/mesh[ixp].dens)*(GAMMA-1.0);
    eneg_R = mesh[ixp].eneg;
    enth_R = (eneg_R + pres_R)/dens_R;
    etrp_R = pres_R/pow(dens_R,GAMMA-1.0);

    REAL ratio = sqrt(dens_R/dens_L);
    REAL vel = (velx_L + ratio*velx_R)/(1.0+ratio);
    REAL enth = (enth_L + ratio*enth_R)/(1.0+ratio);
    REAL cs = sqrt((GAMMA-1.0)*(enth-0.5*SQR(vel)));
    cs_L = sqrt(GAMMA*pres_L/dens_L);
    cs_R = sqrt(GAMMA*pres_R/dens_R);
    
    REAL SL, SR;
    SL = fmin(fmin(velx_L-cs_L, velx_R-cs_R),0.0);
    SR = fmax(fmax(velx_R+cs_R, velx_L+cs_L),0.0);

    if(0.0<SL){
      flux[ix].dens_flux = dens_L*velx_L;
      flux[ix].momx_flux = pres_L + dens_L*SQR(velx_L);
      flux[ix].eneg_flux = (eneg_L+pres_L)*velx_L;
      flux[ix].etrp_flux = etrp_L*velx_L;
    }else if(SL<=0.0 && 0.0<=SR) {
      flux[ix].dens_flux = (SR*momx_L - SL*momx_R + SL*SR*(dens_R-dens_L))/(SR-SL);
      flux[ix].momx_flux = (SR*(pres_L + dens_L*SQR(velx_L)) - SL*(pres_R + dens_R*SQR(velx_R))+SL*SR*(momx_R-momx_L))/(SR-SL);
      flux[ix].eneg_flux = (SR*((eneg_L+pres_L)*velx_L) - SL*((eneg_R+pres_R)*velx_R) +SL*SR*(eneg_R-eneg_L))/(SR-SL);
      flux[ix].etrp_flux = (SR*etrp_L*velx_L - SL*etrp_R*velx_R + SL*SR*(etrp_R-etrp_L))/(SR-SL);
    }else if(SR<0.0) {
      flux[ix].dens_flux = dens_R*velx_R;
      flux[ix].momx_flux = pres_R + dens_R*SQR(velx_R);
      flux[ix].eneg_flux = (eneg_R+pres_R)*velx_R;
      flux[ix].etrp_flux = etrp_R*velx_R;
    }
  }

}

REAL massflux(REAL dens, REAL cs, REAL press, REAL pm)
{
  REAL eps = 1.0e-15;

  REAL gam1 = 0.5*(GAMMA+1.0)/GAMMA;
  REAL gam2 = 0.5*(GAMMA-1.0)/GAMMA;

  REAL flux;

  if(pm/press >= 1.0-eps) {
    flux = dens*cs*sqrt(1.0+gam1*(pm/press-1.0));
  }else{
    flux = dens*cs*gam2*(1.0-pm/press)/(1.0-pow(pm/press,gam2));
  }

  return flux;
}

void sonic(REAL *US1, REAL *US2, REAL *US3, 
	   REAL u1, REAL c1, REAL P1, REAL u2, REAL c2, REAL a1, REAL a2)
{
  REAL R1, R2, us, cs, ps, rs;
  
  R1 = a2/(a2-a1);
  R2 = -a1/(a2-a1);
  us = R1*u1 + R2+u2;
  cs = R1*c1 + R2*c2;
  ps = pow(cs/c1, 2.0*GAMMA/(GAMMA-1.0))*P1;
  rs = GAMMA*ps/SQR(cs);

  *US1 = rs;
  *US2 = rs*us;
  *US3 = ps/(GAMMA-1.0)+0.5*rs*SQR(us);
}


void calc_flux_GODUNOV(struct fluid_1d *mesh, struct fluid_flux *flux, 
		       struct run_param *this_run)
{
  REAL gam = (GAMMA+1.0)/(GAMMA-1.0);
  REAL gam2 = 0.5*(GAMMA-1.0)/GAMMA;

  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL dens_L, dens_R;
    REAL velx_L, velx_R;
    REAL momx_L, momx_R;
    REAL eneg_L, eneg_R;
    REAL pres_L, pres_R;
    REAL enth_L, enth_R;
    REAL cs_L, cs_R;

    int ixp=ix+1; if(ixp>this_run->nmesh-1) ixp=this_run->nmesh-1;

    momx_L = mesh[ix].momx;
    velx_L = mesh[ix].momx/mesh[ix].dens;
    dens_L = mesh[ix].dens;
    pres_L = (GAMMA-1.0)*(mesh[ix].eneg-0.5*SQR(mesh[ix].momx)/mesh[ix].dens);
    eneg_L = mesh[ix].eneg;
    enth_L = (eneg_L + pres_L)/dens_L;
    cs_L = sqrt(GAMMA*pres_L/dens_L);

    momx_R = mesh[ixp].momx;
    velx_R = mesh[ixp].momx/mesh[ixp].dens;
    dens_R = mesh[ixp].dens;
    pres_R = (GAMMA-1.0)*(mesh[ixp].eneg-0.5*SQR(mesh[ixp].momx)/mesh[ixp].dens);
    eneg_R = mesh[ixp].eneg;
    enth_R = (eneg_R + pres_R)/dens_R;
    cs_R = sqrt(GAMMA*pres_R/dens_R);

    // supersonic flow to the RIGHT
    if(velx_L >= cs_L){
      flux[ix].dens_flux = dens_L*velx_L;
      flux[ix].momx_flux = pres_L + dens_L*SQR(velx_L);
      flux[ix].eneg_flux = (eneg_L+pres_L)*velx_L;
    }else if(velx_R <= -cs_R){  // supersonic flow to the LEFT
      flux[ix].dens_flux = dens_R*velx_R;
      flux[ix].momx_flux = pres_R + dens_R*SQR(velx_R);
      flux[ix].eneg_flux = (eneg_R+pres_R)*velx_R;
    }else{

      // initial solution
      REAL pm1 = 0.5*(pres_L + pres_R);

      int k=0;
      int kmax=100;
      REAL tol = 1.0e-5;
      REAL pm2;
      REAL mL, mR, vm;
      do {
	if(k != 0) pm1 = pm2;
	mL = massflux(dens_L, cs_L, pres_L, pm1);
	mR = massflux(dens_R, cs_R, pres_R, pm1);
	pm2 = (mL*pres_R+mR*pres_L-mL*mR*(velx_R-velx_L))/(mL+mR);

	k++;
	if(k>kmax) {
	  printf("fixed-point iteration did not converge\n");
	  exit(EXIT_FAILURE);
	}
      }while(fabs(pm2-pm1)>tol);
      
      mL = massflux(dens_L, cs_L, pres_L, pm2);
      mR = massflux(dens_R, cs_R, pres_R, pm2);
      vm = (mL*velx_L+mR*velx_R-(pres_R-pres_L))/(mL+mR);

      // density in the middle
      REAL densm_L, densm_R;
      // left 
      if(pm2 >= pres_L){
	densm_L = dens_L*(1.0+gam*pm2/pres_L)/(gam+pm2/pres_L);
      }else{
	densm_L = dens_L*pow(pm2/pres_L,1.0/GAMMA);
      }
      // right
      if(pm2 >= pres_R){
	densm_R = dens_R*(1.0+gam*pm2/pres_R)/(gam+pm2/pres_R);
      }else{
	densm_R = dens_R*pow(pm2/pres_R,1.0/GAMMA);
      }

      // contact wave 
      REAL densm_I;
      if (vm > 0.0) {
	densm_I = densm_L;
      }else{
	densm_I = densm_R;
      }
      
      // wave speed at the interface x/t = 0
      REAL csm_L = sqrt(GAMMA*pm2/densm_L);
      REAL csm_R = sqrt(GAMMA*pm2/densm_R);
      REAL Sm_L = vm-csm_L;
      REAL Sm_R = vm+csm_R;

      // sonic case
      REAL Um2, Um3;
      if(Sm_L<=0.0 && Sm_R>=0.0) {
	Um2 = densm_I*vm;
	Um3 = pm2/(GAMMA-1.0)+0.5*densm_I*SQR(vm);
      }else if(Sm_L>0.0 && (velx_L - cs_L) < 0.0){
	sonic(&densm_I, &Um2, &Um3, 
	      velx_L, cs_L, pres_L, vm, csm_L, velx_L-cs_L, Sm_L);
      }else if(Sm_R<0.0 && (velx_R + cs_R) > 0.0){
	sonic(&densm_I, &Um2, &Um3, 
	      velx_R, cs_R, pres_R, vm, csm_R, velx_R+cs_R, Sm_R);
      }

      REAL pres_m = (GAMMA-1.0)*(Um3 - 0.5*SQR(Um2)/densm_I);
      flux[ix].dens_flux = Um2;
      flux[ix].momx_flux = SQR(Um2)/densm_I + pres_m;
      flux[ix].eneg_flux = (Um3 + pres_m)*Um2/densm_I;

    }

  }
}

struct fluid_flux physical_flux(struct fluid_1d *mesh)
{
  struct fluid_flux flux;

  REAL pres = 
    (GAMMA-1.0)*(mesh->eneg-0.5*SQR(mesh->momx)/mesh->dens);
  flux.dens_flux = mesh->momx;
  flux.momx_flux = SQR(mesh->momx)/mesh->dens + pres;
  flux.eneg_flux = (mesh->eneg + pres)*mesh->momx/mesh->dens;
}

struct fluid_flux LF_flux(struct fluid_1d *mi, struct fluid_1d *mip, 
			  struct run_param *this_run)
{
  struct fluid_flux flux, flux_i, flux_ip;
  REAL lambda = this_run->dtime/this_run->delta_x;

  flux_i  = physical_flux(mi);
  flux_ip = physical_flux(mip);

  flux.dens_flux = 0.5*(flux_i.dens_flux + flux_ip.dens_flux + 
			(this_run->CFL/lambda)*(mi->dens - mip->dens));
  flux.momx_flux = 0.5*(flux_i.momx_flux + flux_ip.momx_flux + 
			(this_run->CFL/lambda)*(mi->momx - mip->momx));
  flux.eneg_flux = 0.5*(flux_i.eneg_flux + flux_ip.eneg_flux + 
			(this_run->CFL/lambda)*(mi->eneg - mip->eneg));
  
}

REAL pres_of_mesh(struct fluid_1d *mesh)
{
  REAL pres;

  pres = (GAMMA-1.0)*(mesh->eneg - 0.5*SQR(mesh->momx)/mesh->dens);

  return pres;
}

void positivity_limiter(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
			struct run_param *this_run)
{
  // positivity preserving schemes by
  // Hu, Adams & Shu 2013, Journal of Computational Physics, 242, 169

  REAL eps_dens, eps_pres;
  struct fluid_flux LF_flux_iph;

#if 0
  eps_dens = eps_pres = FLT_MAX;
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL pres 
      = (GAMMA-1.0)*(mesh[ix].eneg-0.5*SQR(mesh[ix].momx)/mesh[ix].dens);
    eps_dens = fmin(eps_dens, mesh[ix].dens);
    eps_pres = fmin(eps_pres, pres);
  }
  assert(eps_dens > 0.0);
  assert(eps_pres > 0.0);
#else
  eps_dens = eps_pres = 1.0e-13;
#endif

  REAL lambda = this_run->dtime/this_run->delta_x;

  for(int ix=0;ix<this_run->nmesh-1;ix++) {
    REAL theta_p, theta_m;

    // positivity for density
    LF_flux_iph = LF_flux(mesh+ix, mesh+ix+1, this_run);

    theta_p = theta_m = 1.0;

    REAL dens_plus  = mesh[ix].dens - 2.0*lambda*flux_iph[ix].dens_flux;
    REAL dens_minus = mesh[ix+1].dens + 2.0*lambda*flux_iph[ix].dens_flux;

    REAL dens_LF_plus  = mesh[ix].dens - 2.0*lambda*LF_flux_iph.dens_flux;
    REAL dens_LF_minus = mesh[ix+1].dens + 2.0*lambda*LF_flux_iph.dens_flux;

    if( dens_plus < eps_dens ) {
      theta_p = (eps_dens - dens_LF_plus)/(dens_plus - dens_LF_plus);
    }

    if( dens_minus < eps_dens) {
      theta_m = (eps_dens - dens_LF_minus)/(dens_minus - dens_LF_minus);
    }

    REAL theta_dens = fmin(theta_p, theta_m);

    flux_iph[ix].dens_flux = 
      (1.0-theta_dens)*LF_flux_iph.dens_flux + 
      theta_dens*flux_iph[ix].dens_flux;

    flux_iph[ix].momx_flux = 
      (1.0-theta_dens)*LF_flux_iph.momx_flux + 
      theta_dens*flux_iph[ix].momx_flux;

    flux_iph[ix].eneg_flux = 
      (1.0-theta_dens)*LF_flux_iph.eneg_flux + 
      theta_dens*flux_iph[ix].eneg_flux;
    

    // positivity for pressure
    theta_p = theta_m = 1.0;

    struct fluid_1d mesh_plus, mesh_minus;
    struct fluid_1d mesh_LF_plus, mesh_LF_minus;

    mesh_plus.dens = mesh[ix].dens - 2.0*lambda*flux_iph[ix].dens_flux;
    mesh_plus.momx = mesh[ix].momx - 2.0*lambda*flux_iph[ix].momx_flux;
    mesh_plus.eneg = mesh[ix].eneg - 2.0*lambda*flux_iph[ix].eneg_flux;

    mesh_minus.dens = mesh[ix+1].dens + 2.0*lambda*flux_iph[ix].dens_flux;
    mesh_minus.momx = mesh[ix+1].momx + 2.0*lambda*flux_iph[ix].momx_flux;
    mesh_minus.eneg = mesh[ix+1].eneg + 2.0*lambda*flux_iph[ix].eneg_flux;

    mesh_LF_plus.dens = mesh[ix].dens - 2.0*lambda*LF_flux_iph.dens_flux;
    mesh_LF_plus.momx = mesh[ix].momx - 2.0*lambda*LF_flux_iph.momx_flux;
    mesh_LF_plus.eneg = mesh[ix].eneg - 2.0*lambda*LF_flux_iph.eneg_flux;

    mesh_LF_minus.dens = mesh[ix+1].dens + 2.0*lambda*LF_flux_iph.dens_flux;
    mesh_LF_minus.momx = mesh[ix+1].momx + 2.0*lambda*LF_flux_iph.momx_flux;
    mesh_LF_minus.eneg = mesh[ix+1].eneg + 2.0*lambda*LF_flux_iph.eneg_flux;

    REAL pres_plus  = pres_of_mesh(&mesh_plus);
    REAL pres_minus = pres_of_mesh(&mesh_minus);
    
    REAL pres_LF_plus  = pres_of_mesh(&mesh_LF_plus);
    REAL pres_LF_minus = pres_of_mesh(&mesh_LF_minus);
    
    if(pres_plus < eps_pres) {
      theta_p = (eps_pres - pres_LF_plus)/(pres_plus - pres_LF_plus);
    }
    if(pres_minus < eps_pres) {
      theta_m = (eps_pres - pres_LF_minus)/(pres_minus - pres_LF_minus);
    }

    REAL theta_pres = fmin(theta_p, theta_m);

    flux_iph[ix].dens_flux = 
      (1.0-theta_pres)*LF_flux_iph.dens_flux + 
      theta_pres*flux_iph[ix].dens_flux;

    flux_iph[ix].momx_flux = 
      (1.0-theta_pres)*LF_flux_iph.momx_flux + 
      theta_pres*flux_iph[ix].momx_flux;

    flux_iph[ix].eneg_flux = 
      (1.0-theta_pres)*LF_flux_iph.eneg_flux + 
      theta_pres*flux_iph[ix].eneg_flux;
  }
  
  
}

void calc_flux(struct fluid_1d *mesh, struct fluid_flux *flux, 
	       struct run_param *this_run)
{
#ifdef __FVS_MUSCL__
  calc_flux_FVS_MUSCL(mesh, flux, this_run);
#elif __FVS_MUSCL2__
  calc_flux_FVS_MUSCL2(mesh, flux, this_run);
#elif __FVS_MP5__
  calc_flux_FVS_MP5(mesh, flux, this_run);
#elif __FVS_MP5_2__
  calc_flux_FVS_MP5_2(mesh, flux, this_run);
#elif __FVS_MP5SL__
  calc_flux_FVS_MP5SL(mesh, flux, this_run);
#elif __AUSMP__
  calc_flux_AUSMP(mesh, flux, this_run);
#elif __AUSMP_MUSCL__
  calc_flux_AUSMP_MUSCL(mesh, flux, this_run);
#elif __AUSMP_MP5__
  calc_flux_AUSMP_MP5(mesh, flux, this_run);
#elif __HLL__
  calc_flux_HLL(mesh, flux, this_run);
#elif __HLLC__
  calc_flux_HLLC(mesh, flux, this_run);
#elif __HLLC_MUSCL__
  calc_flux_HLLC_MUSCL(mesh, flux, this_run);
#elif __HLLC_MP5__
  calc_flux_HLLC_MP5(mesh, flux, this_run);
#elif __GODUNOV__
  calc_flux_GODUNOV(mesh, flux, this_run);
#else
  calc_flux_FVS(mesh, flux, this_run);
#endif

#ifdef __POSITIVITY_LIMITER__
  positivity_limiter(mesh, flux, this_run);
#endif
}
