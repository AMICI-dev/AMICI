
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_w.h"

int dJydp_model_jakstat_adjoint(realtype t, int it, realtype *dJydp, realtype *y, N_Vector x, realtype *dydp, realtype *my, realtype *sigma_y, realtype *dsigma_ydp, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(dJydp,0,sizeof(realtype)*udata->nytrue*udata->nplist*udata->ng);
status = w_model_jakstat_adjoint(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0* udata->nt+it])){
    iy = 0;
  dJydp[iy+(0*17+4)*3] = dydp[12]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+11)*3] = dydp[33]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+13)*3] = dydp[39]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+14)*3] = -dsigma_ydp[42]*(1.0/(sigma_y[0]*sigma_y[0]*sigma_y[0])*pow(my[it+udata->nt*0]-y[it+udata->nt*0],2.0)*1.0-1.0/sigma_y[0]);
}
if(!amiIsNaN(my[1* udata->nt+it])){
    iy = 1;
  dJydp[iy+(0*17+4)*3] = dydp[13]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+10)*3] = dydp[31]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+12)*3] = dydp[37]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+15)*3] = -dsigma_ydp[46]*(1.0/(sigma_y[1]*sigma_y[1]*sigma_y[1])*pow(my[it+udata->nt*1]-y[it+udata->nt*1],2.0)*1.0-1.0/sigma_y[1]);
}
if(!amiIsNaN(my[2* udata->nt+it])){
    iy = 2;
  dJydp[iy+(0*17+5)*3] = dydp[17]*1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+6)*3] = dydp[20]*1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+7)*3] = dydp[23]*1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+8)*3] = dydp[26]*1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+9)*3] = dydp[29]*1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydp[iy+(0*17+16)*3] = -dsigma_ydp[50]*(1.0/(sigma_y[2]*sigma_y[2]*sigma_y[2])*pow(my[it+udata->nt*2]-y[it+udata->nt*2],2.0)*1.0-1.0/sigma_y[2]);
}
return(status);

}


