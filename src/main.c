#include "run_param.h"
#include "fluid.h"

void setup_IC(struct fluid_1d *mesh, struct run_param *this_run)
{

#if __SOD_INIT__
  /* Sod problem */
  sprintf(this_run->model_name,"sod");
  this_run->step = 0;
  this_run->tnow = 0.0;
  this_run->tend = 0.1;
  this_run->xmin = 0.0;
  this_run->xmax = 1.0;
  this_run->delta_x = (this_run->xmax-this_run->xmin)/(REAL)(this_run->nmesh);
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL x = this_run->xmin + ((REAL)ix+0.5)*this_run->delta_x;
    if(x<0.5) {
      REAL press = 1.0;
      mesh[ix].dens = 1.0;
      mesh[ix].momx = 0.0;
      mesh[ix].eneg = press/(GAMMA-1.0);
      mesh[ix].etrp = press/pow(mesh[ix].dens, GAMMA-1.0);
    }else{
      REAL press = 0.1;
      mesh[ix].dens = 0.125;
      mesh[ix].momx = 0.0;
      mesh[ix].eneg = press/(GAMMA-1.0);
      mesh[ix].etrp = press/pow(mesh[ix].dens, GAMMA-1.0);
    }
  }
#elif __LAX_INIT__
  /* Lax problem */
  sprintf(this_run->model_name,"lax");
  this_run->step = 0;
  this_run->tnow = 0.0;
  this_run->tend = 0.32;
  this_run->xmin =-1.0;
  this_run->xmax = 1.0;
  this_run->delta_x = (this_run->xmax-this_run->xmin)/(REAL)(this_run->nmesh);
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL x = this_run->xmin + ((REAL)ix+0.5)*this_run->delta_x;
    if(x<0.0) {
      REAL press = 3.527;
      mesh[ix].dens = 0.445;
      mesh[ix].momx = 0.698*mesh[ix].dens;
      mesh[ix].eneg = 0.5*SQR(mesh[ix].momx)/mesh[ix].dens + press/(GAMMA-1.0);
      mesh[ix].etrp = press/pow(mesh[ix].dens, GAMMA-1.0);
    }else{
      REAL press = 0.571;
      mesh[ix].dens = 0.5;
      mesh[ix].momx = 0.0;
      mesh[ix].eneg = 0.5*SQR(mesh[ix].momx)/mesh[ix].dens + press/(GAMMA-1.0);
      mesh[ix].etrp = press/pow(mesh[ix].dens, GAMMA-1.0);
    }
  }
#elif __SHU_INIT__
  /* Shu's problem */
  sprintf(this_run->model_name,"shu");
  this_run->step = 0;
  this_run->tnow = 0.0;
  this_run->tend = 0.36;
  this_run->xmin =-1.0;
  this_run->xmax = 1.0;
  this_run->delta_x = (this_run->xmax-this_run->xmin)/(REAL)(this_run->nmesh);
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL x = this_run->xmin + ((REAL)ix+0.5)*this_run->delta_x;
    if(x<-0.8) {
      REAL press = 10.3333;
      mesh[ix].dens = 3.857143;
      mesh[ix].momx = 2.629369*mesh[ix].dens;
      mesh[ix].eneg = 0.5*SQR(mesh[ix].momx)/mesh[ix].dens + press/(GAMMA-1.0);
      mesh[ix].etrp = press/pow(mesh[ix].dens, GAMMA-1.0);
    }else{
      REAL press = 1.0;
      mesh[ix].dens = 1.0+0.2*sin(5.0*PI*x);
      mesh[ix].momx = 0.0;
      mesh[ix].eneg = 0.5*SQR(mesh[ix].momx)/mesh[ix].dens + press/(GAMMA-1.0);
      mesh[ix].etrp = press/pow(mesh[ix].dens, GAMMA-1.0);
    }

  }
#elif __LE_BLANC_INIT__
  /* Le Blanc problem */
  sprintf(this_run->model_name,"le_blanc");
  this_run->step = 0;
  this_run->tnow = 0.0;
  this_run->tend = 6.0;
  this_run->xmin = 0.0;
  this_run->xmax = 9.0;
  this_run->delta_x = (this_run->xmax-this_run->xmin)/(REAL)(this_run->nmesh);
  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL x = this_run->xmin + ((REAL)ix+0.5)*this_run->delta_x;
    if(x<3.0) {
      REAL press = 0.2/3.0;
      mesh[ix].dens = 1.0;
      mesh[ix].momx = 0.0;
      mesh[ix].eneg = 0.5*SQR(mesh[ix].momx)/mesh[ix].dens + press/(GAMMA-1.0);
      mesh[ix].etrp = press/pow(mesh[ix].dens, GAMMA-1.0);
    }else{
      REAL press = 2.0e-10/3.0;
      mesh[ix].dens = 1.0e-3;
      mesh[ix].momx = 0.0;
      mesh[ix].eneg = 0.5*SQR(mesh[ix].momx)/mesh[ix].dens + press/(GAMMA-1.0);
      mesh[ix].etrp = press/pow(mesh[ix].dens, GAMMA-1.0);      
    }

  }
#else
#error "NO INIT FLAG SPECIFIED"
#endif

#ifdef __FVS_MUSCL__
  sprintf(this_run->scheme_name,"FVS_MUSCL");
#elif __FVS_MUSCL2__
  sprintf(this_run->scheme_name,"FVS_MUSCL2");
#elif __FVS_MP5__
#ifdef __TVD_RK3__
  sprintf(this_run->scheme_name,"FVS_MP5_TVDRK3");
#else
  sprintf(this_run->scheme_name,"FVS_MP5_WO_TVDRK3");
#endif
#elif __FVS_MP5_2__
#ifdef __TVD_RK3__
  sprintf(this_run->scheme_name,"FVS_MP5_2_TVDRK3");
#else
  sprintf(this_run->scheme_name,"FVS_MP5_2_WO_TVDRK3");
#endif
#elif __FVS_MP5SL__
  sprintf(this_run->scheme_name,"FVS_MP5SL");
#elif __AUSMP__
  sprintf(this_run->scheme_name,"AUSMP");
#elif __AUSMP_MUSCL__
  sprintf(this_run->scheme_name,"AUSMP_MUSCL");
#elif __AUSMP_MP5__
  sprintf(this_run->scheme_name,"AUSMP_MP5");
#elif __HLL__
  sprintf(this_run->scheme_name,"HLL");
#elif __HLLC__
  sprintf(this_run->scheme_name,"HLLC");
#elif __HLLC_MUSCL__
  sprintf(this_run->scheme_name,"HLLC_MUSCL");
#elif __HLLC_MP5__
  sprintf(this_run->scheme_name,"HLLC_MP5");
#else
  sprintf(this_run->scheme_name,"FVS");
#endif
}

void sweep_1d(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
	      struct run_param *this_run)
{
  REAL nu = this_run->dtime/this_run->delta_x;

#ifdef __PERIODIC__
  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ixm = (ix-1+this_run->nmesh)%this_run->nmesh;
    mesh[ix].dens -= nu*(flux_iph[ix].dens_flux - flux_iph[ixm].dens_flux);
    mesh[ix].momx -= nu*(flux_iph[ix].momx_flux - flux_iph[ixm].momx_flux);
    mesh[ix].eneg -= nu*(flux_iph[ix].eneg_flux - flux_iph[ixm].eneg_flux);
#if defined(__ENTROPY_INTEGRATION_AVAILABLE__) && defined(__ENTROPY_INTEGRATION__)
    REAL eneg_k = 0.5*SQR(mesh[ix].momx)/mesh[ix].dens;
    if(mesh[ix].high_mach==1 && mesh[ix].under_shock==0) {
      mesh[ix].etrp -= nu*(flux_iph[ix].etrp_flux - flux_iph[ixm].etrp_flux);
      REAL eneg_th = mesh[ix].etrp*pow(mesh[ix].dens,GAMMA-1.0)/(GAMMA-1.0);
      mesh[ix].eneg =  eneg_k + eneg_th;
    }else{
      REAL eneg_th = mesh[ix].eneg - eneg_k;
      mesh[ix].etrp = (GAMMA-1.0)*eneg_th/pow(mesh[ix].dens, GAMMA-1.0);
    }
#endif
  }  
#else
  for(int ix=1;ix<this_run->nmesh-1;ix++) {
    int ixm = ix-1;
    mesh[ix].dens -= nu*(flux_iph[ix].dens_flux - flux_iph[ixm].dens_flux);
    mesh[ix].momx -= nu*(flux_iph[ix].momx_flux - flux_iph[ixm].momx_flux);
    mesh[ix].eneg -= nu*(flux_iph[ix].eneg_flux - flux_iph[ixm].eneg_flux);
#if defined(__ENTROPY_INTEGRATION_AVAILABLE__) && defined(__ENTROPY_INTEGRATION__)
    REAL eneg_k = 0.5*SQR(mesh[ix].momx)/mesh[ix].dens;
    if(mesh[ix].high_mach==1 && mesh[ix].under_shock==0) {
      mesh[ix].etrp -= nu*(flux_iph[ix].etrp_flux - flux_iph[ixm].etrp_flux);
      REAL eneg_th = mesh[ix].etrp*pow(mesh[ix].dens,GAMMA-1.0)/(GAMMA-1.0);
      mesh[ix].eneg =  eneg_k + eneg_th;
    }else{
      REAL eneg_th = mesh[ix].eneg - eneg_k;
      mesh[ix].etrp = (GAMMA-1.0)*eneg_th/pow(mesh[ix].dens, GAMMA-1.0);
    }
#endif
  }
#endif
}

void integrate(struct fluid_1d *mesh, struct fluid_flux *flux_iph, 
	       struct run_param *this_run)
{
#ifdef __TVD_RK3__
  struct fluid_1d *mesh_old;

  mesh_old = (struct fluid_1d*) malloc(sizeof(struct fluid_1d)*this_run->nmesh);
  for(int ix=0;ix<this_run->nmesh;ix++) mesh_old[ix] = mesh[ix];

  calc_flux(mesh, flux_iph, this_run);
  sweep_1d(mesh, flux_iph, this_run);

  calc_flux(mesh, flux_iph, this_run);
  sweep_1d(mesh, flux_iph, this_run);
  for(int ix=0;ix<this_run->nmesh;ix++) {
    mesh[ix].dens = 0.75*mesh_old[ix].dens + 0.25*mesh[ix].dens;
    mesh[ix].momx = 0.75*mesh_old[ix].momx + 0.25*mesh[ix].momx;
    mesh[ix].eneg = 0.75*mesh_old[ix].eneg + 0.25*mesh[ix].eneg;
  }

  calc_flux(mesh, flux_iph, this_run);
  sweep_1d(mesh, flux_iph, this_run);
  for(int ix=0;ix<this_run->nmesh;ix++) {
    mesh[ix].dens = 1.0/3.0*mesh_old[ix].dens + 2.0/3.0*mesh[ix].dens;
    mesh[ix].momx = 1.0/3.0*mesh_old[ix].momx + 2.0/3.0*mesh[ix].momx;
    mesh[ix].eneg = 1.0/3.0*mesh_old[ix].eneg + 2.0/3.0*mesh[ix].eneg;
  }

  free(mesh_old);
#else /* !__TVD_RK3__ */
  calc_flux(mesh, flux_iph, this_run);
  sweep_1d(mesh, flux_iph, this_run);
#endif /* __TVD_RK3__ */

  this_run->tnow += this_run->dtime;
  this_run->step++;
}

void calc_dtime(struct fluid_1d *mesh, struct run_param *this_run)
{
  this_run->dtime = FLT_MAX;

  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL eneg_th = mesh[ix].eneg - 0.5*SQR(mesh[ix].momx)/mesh[ix].dens;
    REAL cs = sqrt(GAMMA*(GAMMA-1.0)*eneg_th/mesh[ix].dens);
    REAL velx = fabs(mesh[ix].momx/mesh[ix].dens);

    this_run->dtime = fmin(this_run->dtime, this_run->delta_x/(velx+cs));
  }
  this_run->dtime *= this_run->CFL;

  if(this_run->tnow + this_run->dtime > this_run->tend) this_run->dtime = this_run->tend - this_run->tnow;

}

void output_fluid(struct fluid_1d *mesh, struct run_param *this_run, FILE *fp)
{

  for(int ix=0;ix<this_run->nmesh;ix++) {
    REAL x = this_run->xmin + ((REAL)ix+0.5)*this_run->delta_x;
    REAL ene_th = mesh[ix].eneg - 0.5*SQR(mesh[ix].momx)/mesh[ix].dens;
    fprintf(fp, "%d %14.6e %14.6e %14.6e %14.6e %d %d\n", ix, x, 
	    mesh[ix].dens, mesh[ix].momx/mesh[ix].dens, ene_th*(GAMMA-1.0),
	    mesh[ix].high_mach, mesh[ix].under_shock);
  }
}

void output_fluid_with_name(struct fluid_1d *mesh, struct run_param *this_run)
{
  FILE *fp;

  static char filename[64];
  
  sprintf(filename, "%s_%s_%05d.dat", this_run->model_name, this_run->scheme_name, this_run->step);
  fp = fopen(filename, "w");
  output_fluid(mesh, this_run, fp);
  fclose(fp);
}

void input_fluid(struct fluid_1d *mesh, struct run_param *this_run,
		 char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");

  for(int ix=0;ix<this_run->nmesh;ix++) {
    int ix;
    REAL x, dens, velx, pres;
    fscanf(fp, "%d %f %f %f %f", &ix, &x, &dens, &velx, &pres);
    mesh[ix].dens = dens;
    mesh[ix].momx = dens*velx;
    mesh[ix].eneg = 0.5*dens*SQR(velx) + pres/(GAMMA-1.0);
  }

  fclose(fp);
}

int main(int argc, char **argv)
{
  struct run_param this_run;
  struct fluid_1d *mesh;
  struct fluid_flux *flux;

  this_run.nmesh = __NMESH__;
  this_run.CFL = 0.2;
  
  mesh = (struct fluid_1d *) malloc(sizeof(struct fluid_1d)*this_run.nmesh);
  flux = (struct fluid_flux *) malloc(sizeof(struct fluid_flux)*this_run.nmesh);
  
  setup_IC(mesh, &this_run);

  calc_dtime(mesh, &this_run);

  while(this_run.tnow < this_run.tend) {
    integrate(mesh, flux, &this_run);
    calc_dtime(mesh, &this_run); 
    //    if((this_run.step%10) == 0) output_fluid_with_name(mesh, &this_run);
  }

  output_fluid(mesh, &this_run, stdout);
  
}
