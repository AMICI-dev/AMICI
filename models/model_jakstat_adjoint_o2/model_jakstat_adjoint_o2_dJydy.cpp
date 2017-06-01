
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

int dJydy_model_jakstat_adjoint_o2(realtype t, int it, realtype *dJydy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(dJydy,0,sizeof(realtype)*udata->nytrue*udata->nytrue*udata->ng);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0* udata->nt+it])){
    iy = 0;
  dJydy[iy+(0*18+0)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+1)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*3]*1.0;
  dJydy[iy+(3*18+1)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+2)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*6]*1.0;
  dJydy[iy+(6*18+2)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+3)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*9]*1.0;
  dJydy[iy+(9*18+3)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+4)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*12]*1.0;
  dJydy[iy+(12*18+4)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+5)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*15]*1.0;
  dJydy[iy+(15*18+5)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+6)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*18]*1.0;
  dJydy[iy+(18*18+6)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+7)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*21]*1.0;
  dJydy[iy+(21*18+7)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+8)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*24]*1.0;
  dJydy[iy+(24*18+8)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+9)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*27]*1.0;
  dJydy[iy+(27*18+9)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+10)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*30]*1.0;
  dJydy[iy+(30*18+10)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+11)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*33]*1.0;
  dJydy[iy+(33*18+11)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+12)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*36]*1.0;
  dJydy[iy+(36*18+12)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+13)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*39]*1.0;
  dJydy[iy+(39*18+13)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+14)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*42]*1.0;
  dJydy[iy+(42*18+14)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+15)*3] = 1.0/(sigma_y[0]*sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*1.0+1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*45]*1.0;
  dJydy[iy+(45*18+15)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+16)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*48]*1.0;
  dJydy[iy+(48*18+16)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  dJydy[iy+(0*18+17)*3] = 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*51]*1.0;
  dJydy[iy+(51*18+17)*3] = 1.0/(sigma_y[0]*sigma_y[0])*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[1* udata->nt+it])){
    iy = 1;
  dJydy[iy+(1*18+0)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+1)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*4]*1.0;
  dJydy[iy+(4*18+1)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+2)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*7]*1.0;
  dJydy[iy+(7*18+2)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+3)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*10]*1.0;
  dJydy[iy+(10*18+3)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+4)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*13]*1.0;
  dJydy[iy+(13*18+4)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+5)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*16]*1.0;
  dJydy[iy+(16*18+5)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+6)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*19]*1.0;
  dJydy[iy+(19*18+6)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+7)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*22]*1.0;
  dJydy[iy+(22*18+7)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+8)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*25]*1.0;
  dJydy[iy+(25*18+8)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+9)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*28]*1.0;
  dJydy[iy+(28*18+9)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+10)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*31]*1.0;
  dJydy[iy+(31*18+10)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+11)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*34]*1.0;
  dJydy[iy+(34*18+11)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+12)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*37]*1.0;
  dJydy[iy+(37*18+12)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+13)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*40]*1.0;
  dJydy[iy+(40*18+13)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+14)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*43]*1.0;
  dJydy[iy+(43*18+14)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+15)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*46]*1.0;
  dJydy[iy+(46*18+15)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+16)*3] = 1.0/(sigma_y[1]*sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*1.0+1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*49]*1.0;
  dJydy[iy+(49*18+16)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  dJydy[iy+(1*18+17)*3] = 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*52]*1.0;
  dJydy[iy+(52*18+17)*3] = 1.0/(sigma_y[1]*sigma_y[1])*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[2* udata->nt+it])){
    iy = 2;
  dJydy[iy+(2*18+0)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+1)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*5]*1.0;
  dJydy[iy+(5*18+1)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+2)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*8]*1.0;
  dJydy[iy+(8*18+2)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+3)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*11]*1.0;
  dJydy[iy+(11*18+3)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+4)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*14]*1.0;
  dJydy[iy+(14*18+4)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+5)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*17]*1.0;
  dJydy[iy+(17*18+5)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+6)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*20]*1.0;
  dJydy[iy+(20*18+6)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+7)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*23]*1.0;
  dJydy[iy+(23*18+7)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+8)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*26]*1.0;
  dJydy[iy+(26*18+8)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+9)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*29]*1.0;
  dJydy[iy+(29*18+9)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+10)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*32]*1.0;
  dJydy[iy+(32*18+10)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+11)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*35]*1.0;
  dJydy[iy+(35*18+11)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+12)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*38]*1.0;
  dJydy[iy+(38*18+12)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+13)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*41]*1.0;
  dJydy[iy+(41*18+13)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+14)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*44]*1.0;
  dJydy[iy+(44*18+14)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+15)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*47]*1.0;
  dJydy[iy+(47*18+15)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+16)*3] = 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*50]*1.0;
  dJydy[iy+(50*18+16)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  dJydy[iy+(2*18+17)*3] = 1.0/(sigma_y[2]*sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*1.0+1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*53]*1.0;
  dJydy[iy+(53*18+17)*3] = 1.0/(sigma_y[2]*sigma_y[2])*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
}
return(status);

}


