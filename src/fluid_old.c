#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define NMESH (128)
#define GAMMA (1.4)
#define COURANT (0.2)

#define SQR(a) ((a)*(a))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

struct fluid1D{
  float rho;
  float mom;
  float ene;
};

void setup_IC(struct fluid1D *fluid)
{
  int i;

  for(i=1;i<=NMESH;i++) {
    float dens;
    float pres;
    float uene;

    if(i<NMESH/2) {
      dens = 1.0;
      pres = 1.0;
    }else{
      dens = 0.125;
      pres = 0.1;
    }
    uene = pres/(GAMMA-1.0);

    fluid[i].rho=dens;
    fluid[i].mom=0.0;
    fluid[i].ene=0.5*SQR(fluid[i].mom)/fluid[i].rho+uene;
  }
  
  return;
}

void sweep_bc(struct fluid1D* fluid)
{
  fluid[0]=fluid[1];
  fluid[NMESH+1]=fluid[NMESH];

  return;
}

void sweep_1d(struct fluid1D *fluid, float dt)
{
  int i;
  int k,l,m;

  float dx;

  static  float Fp[NMESH+2][3];
  static  float Fm[NMESH+2][3];

  static float lambda[3];
  static float Rm[9], Lm[9], Wpm[9], Wmm[9];

#define R(a,b) Rm[(a)+3*(b)] /* 0 <= a,b <= 2 */
#define L(a,b) Lm[(a)+3*(b)] /* 0 <= a,b <= 2 */
#define Wp(a,b) Wpm[(a)+3*(b)] /* 0 <= a,b <= 2 */
#define Wm(a,b) Wmm[(a)+3*(b)] /* 0 <= a,b <= 2 */

  dx = 1.0/(float)NMESH;
  
  for(i=0;i<NMESH+2;i++) {
    float uene, pres, vel, enth, cs;
    float b1, b2;

    vel  = fluid[i].mom/fluid[i].rho;
    uene = fluid[i].ene-0.5*fluid[i].rho*SQR(vel);
    pres = uene*(GAMMA-1.0);
    enth = (fluid[i].ene+pres)/fluid[i].rho;
    cs = sqrt(GAMMA*pres/fluid[i].rho);

    lambda[0]=vel-cs;
    lambda[1]=vel;
    lambda[2]=vel+cs;
    
    R(0,0) = 1.0;
    R(0,1) = 1.0;
    R(0,2) = 1.0;

    R(1,0) = vel-cs;
    R(1,1) = vel;
    R(1,2) = vel+cs;

    R(2,0) = enth-vel*cs;
    R(2,1) = 0.5*SQR(vel);
    R(2,2) = enth+vel*cs;

    b1 = 0.5*(GAMMA-1.0)*SQR(vel/cs);
    b2 = (GAMMA-1.0)/SQR(cs);

    L(0,0) = 0.5*(b1+vel/cs);
    L(0,1) = -0.5*(1.0/cs+b2*vel);
    L(0,2) = 0.5*b2;

    L(1,0) = 1.0-b1;
    L(1,1) = b2*vel;
    L(1,2) = -b2;

    L(2,0) = 0.5*(b1-vel/cs);
    L(2,1) = 0.5*(1.0/cs-b2*vel);
    L(2,2) = 0.5*b2;

    for(k=0;k<3;k++) {
      for(l=0;l<3;l++) {
	Wp(k,l)=0.0;
	Wm(k,l)=0.0;
	for(m=0;m<3;m++) Wp(k,l) += R(k,m)*MAX(lambda[m],0.0)*L(m,l);
	for(m=0;m<3;m++) Wm(k,l) += R(k,m)*MIN(lambda[m],0.0)*L(m,l);
      }
    }

    Fp[i][0] = Wp(0,0)*fluid[i].rho+Wp(0,1)*fluid[i].mom+Wp(0,2)*fluid[i].ene;
    Fp[i][1] = Wp(1,0)*fluid[i].rho+Wp(1,1)*fluid[i].mom+Wp(1,2)*fluid[i].ene;
    Fp[i][2] = Wp(2,0)*fluid[i].rho+Wp(2,1)*fluid[i].mom+Wp(2,2)*fluid[i].ene;

    Fm[i][0] = Wm(0,0)*fluid[i].rho+Wm(0,1)*fluid[i].mom+Wm(0,2)*fluid[i].ene;
    Fm[i][1] = Wm(1,0)*fluid[i].rho+Wm(1,1)*fluid[i].mom+Wm(1,2)*fluid[i].ene;
    Fm[i][2] = Wm(2,0)*fluid[i].rho+Wm(2,1)*fluid[i].mom+Wm(2,2)*fluid[i].ene;

  }

#undef Wm
#undef Wp
#undef L
#undef R

  for(i=1;i<=NMESH;i++) {
    fluid[i].rho = fluid[i].rho 
      - dt/dx*(Fp[i][0]+Fm[i+1][0]-Fp[i-1][0]-Fm[i][0]);
    fluid[i].mom = fluid[i].mom
      - dt/dx*(Fp[i][1]+Fm[i+1][1]-Fp[i-1][1]-Fm[i][1]);
    fluid[i].ene = fluid[i].ene
      - dt/dx*(Fp[i][2]+Fm[i+1][2]-Fp[i-1][2]-Fm[i][2]);
  }
}

#define KAPPA (0.0)

float muscl_L(float u_im1, float u_i, float u_ip1)
{
  float delta_p, delta_m, delta_pp, delta_mm, beta;

  beta = (3.0-KAPPA)/(1.0-KAPPA);

  delta_p = u_ip1-u_i;
  delta_m = u_i-u_im1;

  if(delta_p*beta*delta_m<0.0) {
    delta_pp=0.0;
  }else if(fabs(delta_p) < fabs(beta*delta_m)) {
    delta_pp=delta_p;
  }else{
    delta_pp = beta*delta_m;
  }

  if(delta_m*beta*delta_p<0.0) {
    delta_mm=0.0;
  }else if(fabs(delta_m) < fabs(beta*delta_p)){
    delta_mm=delta_m;
  }else{
    delta_mm=beta*delta_p;
  }

  return (u_i+0.25*(1.0-KAPPA)*delta_mm+0.25*(1.0+KAPPA)*delta_pp);
}

float muscl_R(float u_im1, float u_i, float u_ip1)
{
  float delta_p, delta_m, delta_pp, delta_mm, beta;

  beta = (3.0-KAPPA)/(1.0-KAPPA);

  delta_p = u_ip1-u_i;
  delta_m = u_i-u_im1;

  if(delta_p*beta*delta_m<0.0) {
    delta_pp=0.0;
  }else if(fabs(delta_p) < fabs(beta*delta_m)) {
    delta_pp=delta_p;
  }else{
    delta_pp = beta*delta_m;
  }

  if(delta_m*beta*delta_p<0.0) {
    delta_mm=0.0;
  }else if(fabs(delta_m) < fabs(beta*delta_p)){
    delta_mm=delta_m;
  }else{
    delta_mm=beta*delta_p;
  }

  return (u_i-0.25*(1.0+KAPPA)*delta_mm-0.25*(1.0-KAPPA)*delta_pp);
}

void sweep_1d_muscl(struct fluid1D *fluid, float dt)
{
  int i;
  int k,l,m;

  float dx;

  static  float Fp[NMESH+2][3];
  static  float Fm[NMESH+2][3];

  static float lambda[3];
  static float Rm[9], Lm[9], Wpm[9], Wmm[9];

#define R(a,b) Rm[(a)+3*(b)] /* 0 <= a,b <= 2 */
#define L(a,b) Lm[(a)+3*(b)] /* 0 <= a,b <= 2 */
#define Wp(a,b) Wpm[(a)+3*(b)] /* 0 <= a,b <= 2 */
#define Wm(a,b) Wmm[(a)+3*(b)] /* 0 <= a,b <= 2 */

  dx = 1.0/(float)NMESH;
  
  for(i=0;i<NMESH+2;i++) {
    float rho_L, mom_L, ene_L;
    float rho_R, mom_R, ene_R;

    float uene_L, pres_L, vel_L, enth_L, cs_L;
    float uene_R, pres_R, vel_R, enth_R, cs_R;

    float b1, b2;

    if(i==0) {
      rho_L = muscl_L(fluid[i].rho, fluid[i].rho, fluid[i+1].rho);
      mom_L = muscl_L(fluid[i].mom, fluid[i].mom, fluid[i+1].mom);
      ene_L = muscl_L(fluid[i].ene, fluid[i].ene, fluid[i+1].ene);
      rho_R = muscl_R(fluid[i].rho, fluid[i].rho, fluid[i+1].rho);
      mom_R = muscl_R(fluid[i].mom, fluid[i].mom, fluid[i+1].mom);
      ene_R = muscl_R(fluid[i].ene, fluid[i].ene, fluid[i+1].ene);
    }else if(i==NMESH+1){
      rho_L = muscl_L(fluid[i-1].rho, fluid[i].rho, fluid[i].rho);
      mom_L = muscl_L(fluid[i-1].mom, fluid[i].mom, fluid[i].mom);
      ene_L = muscl_L(fluid[i-1].ene, fluid[i].ene, fluid[i].ene);
      rho_R = muscl_R(fluid[i-1].rho, fluid[i].rho, fluid[i].rho);
      mom_R = muscl_R(fluid[i-1].mom, fluid[i].mom, fluid[i].mom);
      ene_R = muscl_R(fluid[i-1].ene, fluid[i].ene, fluid[i].ene);
    }else{
      rho_L = muscl_L(fluid[i-1].rho, fluid[i].rho, fluid[i+1].rho);
      mom_L = muscl_L(fluid[i-1].mom, fluid[i].mom, fluid[i+1].mom);
      ene_L = muscl_L(fluid[i-1].ene, fluid[i].ene, fluid[i+1].ene);
      rho_R = muscl_R(fluid[i-1].rho, fluid[i].rho, fluid[i+1].rho);
      mom_R = muscl_R(fluid[i-1].mom, fluid[i].mom, fluid[i+1].mom);
      ene_R = muscl_R(fluid[i-1].ene, fluid[i].ene, fluid[i+1].ene);
    }

    vel_L = mom_L/rho_L;
    uene_L = ene_L-0.5*rho_L*SQR(vel_L);
    pres_L = uene_L*(GAMMA-1.0);
    enth_L = (ene_L+pres_L)/rho_L;
    cs_L = sqrt(GAMMA*pres_L/rho_L);

    lambda[0]=vel_L-cs_L;
    lambda[1]=vel_L;
    lambda[2]=vel_L+cs_L;
    
    R(0,0) = 1.0;
    R(0,1) = 1.0;
    R(0,2) = 1.0;

    R(1,0) = vel_L-cs_L;
    R(1,1) = vel_L;
    R(1,2) = vel_L+cs_L;

    R(2,0) = enth_L-vel_L*cs_L;
    R(2,1) = 0.5*SQR(vel_L);
    R(2,2) = enth_L+vel_L*cs_L;

    b1 = 0.5*(GAMMA-1.0)*SQR(vel_L/cs_L);
    b2 = (GAMMA-1.0)/SQR(cs_L);

    L(0,0) = 0.5*(b1+vel_L/cs_L);
    L(0,1) = -0.5*(1.0/cs_L+b2*vel_L);
    L(0,2) = 0.5*b2;

    L(1,0) = 1.0-b1;
    L(1,1) = b2*vel_L;
    L(1,2) = -b2;

    L(2,0) = 0.5*(b1-vel_L/cs_L);
    L(2,1) = 0.5*(1.0/cs_L-b2*vel_L);
    L(2,2) = 0.5*b2;

    for(k=0;k<3;k++) {
      for(l=0;l<3;l++) {
	Wp(k,l)=0.0;
	for(m=0;m<3;m++) Wp(k,l) += R(k,m)*MAX(lambda[m],0.0)*L(m,l);
      }
    }

    vel_R = mom_R/rho_R;
    uene_R = ene_R-0.5*rho_R*SQR(vel_R);
    pres_R = uene_R*(GAMMA-1.0);
    enth_R = (ene_R+pres_R)/rho_R;
    cs_R = sqrt(GAMMA*pres_R/rho_R);

    lambda[0]=vel_R-cs_R;
    lambda[1]=vel_R;
    lambda[2]=vel_R+cs_R;
    
    R(0,0) = 1.0;
    R(0,1) = 1.0;
    R(0,2) = 1.0;

    R(1,0) = vel_R-cs_R;
    R(1,1) = vel_R;
    R(1,2) = vel_R+cs_R;

    R(2,0) = enth_R-vel_R*cs_R;
    R(2,1) = 0.5*SQR(vel_R);
    R(2,2) = enth_R+vel_R*cs_R;

    b1 = 0.5*(GAMMA-1.0)*SQR(vel_R/cs_R);
    b2 = (GAMMA-1.0)/SQR(cs_R);

    L(0,0) = 0.5*(b1+vel_R/cs_R);
    L(0,1) = -0.5*(1.0/cs_R+b2*vel_R);
    L(0,2) = 0.5*b2;

    L(1,0) = 1.0-b1;
    L(1,1) = b2*vel_R;
    L(1,2) = -b2;

    L(2,0) = 0.5*(b1-vel_R/cs_R);
    L(2,1) = 0.5*(1.0/cs_R-b2*vel_R);
    L(2,2) = 0.5*b2;

    for(k=0;k<3;k++) {
      for(l=0;l<3;l++) {
	Wm(k,l)=0.0;
	for(m=0;m<3;m++) Wm(k,l) += R(k,m)*MIN(lambda[m],0.0)*L(m,l);
      }
    }

    Fp[i][0] = Wp(0,0)*rho_L+Wp(0,1)*mom_L+Wp(0,2)*ene_L;
    Fp[i][1] = Wp(1,0)*rho_L+Wp(1,1)*mom_L+Wp(1,2)*ene_L;
    Fp[i][2] = Wp(2,0)*rho_L+Wp(2,1)*mom_L+Wp(2,2)*ene_L;

    Fm[i][0] = Wm(0,0)*rho_R+Wm(0,1)*mom_R+Wm(0,2)*ene_R;
    Fm[i][1] = Wm(1,0)*rho_R+Wm(1,1)*mom_R+Wm(1,2)*ene_R;
    Fm[i][2] = Wm(2,0)*rho_R+Wm(2,1)*mom_R+Wm(2,2)*ene_R;

  }

#undef Wm
#undef Wp
#undef L
#undef R

  for(i=1;i<=NMESH;i++) {
    fluid[i].rho = fluid[i].rho 
      - dt/dx*(Fp[i][0]+Fm[i+1][0]-Fp[i-1][0]-Fm[i][0]);
    fluid[i].mom = fluid[i].mom
      - dt/dx*(Fp[i][1]+Fm[i+1][1]-Fp[i-1][1]-Fm[i][1]);
    fluid[i].ene = fluid[i].ene
      - dt/dx*(Fp[i][2]+Fm[i+1][2]-Fp[i-1][2]-Fm[i][2]);
  }
}

void set_timestep(struct fluid1D *fluid, float *dt)
{
  int i;
  float dtime,dx,vel,uene,pres,cs;

  dtime = FLT_MAX;
  dx = 1.0/(float)NMESH;

  for(i=0;i<NMESH;i++) {
    vel = fluid[i].mom/fluid[i].rho;
    uene = fluid[i].ene-0.5*fluid[i].rho*SQR(vel);
    pres = uene*(GAMMA-1.0);
    cs = sqrt(GAMMA*pres/fluid[i].rho);
    dtime = MIN(dtime, dx/(fabs(vel)+cs));
  }

  *dt = COURANT*dtime;
  
}

int main(int argc, char **argv) 
{

  int i;

  static struct fluid1D fluid[NMESH+2];
  float dt,time;

  setup_IC(fluid);
  sweep_bc(fluid);

  time=0.0;

  while(time<0.1){
    set_timestep(fluid,&dt);
    //    sweep_1d_muscl(fluid, dt);
    sweep_1d(fluid,dt);
    time += dt;
  }

  for(i=1;i<NMESH+1;i++) {
    float vel, uene, pres;
    
    vel = fluid[i].mom/fluid[i].rho;
    uene = fluid[i].ene-0.5*fluid[i].rho*SQR(vel);
    pres = uene*(GAMMA-1.0);
    printf("%14.6e %14.6e %14.6e %14.6e\n", 
	   (float)i/NMESH, fluid[i].rho, vel, pres);
  }

}
