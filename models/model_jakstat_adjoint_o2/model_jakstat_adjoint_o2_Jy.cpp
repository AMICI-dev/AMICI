
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_jakstat_adjoint_o2_w.h"

int Jy_model_jakstat_adjoint_o2(realtype t, int it, N_Vector x, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->Jy[0] += amilog((tdata->sigmay[0]*tdata->sigmay[0])*3.141592653589793*2.0)*5.0E-1+1.0/(tdata->sigmay[0]*tdata->sigmay[0])*pow(edata->my[it+udata->nt*0]-rdata->y[it + udata->nt*0],2.0)*5.0E-1;
  tdata->Jy[1] += rdata->y[it + udata->nt*3]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[2] += rdata->y[it + udata->nt*6]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[3] += rdata->y[it + udata->nt*9]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[4] += rdata->y[it + udata->nt*12]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[5] += rdata->y[it + udata->nt*15]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[6] += rdata->y[it + udata->nt*18]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[7] += rdata->y[it + udata->nt*21]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[8] += rdata->y[it + udata->nt*24]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[9] += rdata->y[it + udata->nt*27]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[10] += rdata->y[it + udata->nt*30]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[11] += rdata->y[it + udata->nt*33]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[12] += rdata->y[it + udata->nt*36]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[13] += rdata->y[it + udata->nt*39]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[14] += rdata->y[it + udata->nt*42]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[15] += 1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*pow(edata->my[it+udata->nt*0]-rdata->y[it + udata->nt*0],2.0)*-1.0+1.0/tdata->sigmay[0]-rdata->y[it + udata->nt*45]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*5.0E-1;
  tdata->Jy[16] += rdata->y[it + udata->nt*48]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->Jy[17] += rdata->y[it + udata->nt*51]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(edata->my[1* udata->nt+it])){
    iy = 1;
  tdata->Jy[0] += amilog((tdata->sigmay[1]*tdata->sigmay[1])*3.141592653589793*2.0)*5.0E-1+1.0/(tdata->sigmay[1]*tdata->sigmay[1])*pow(edata->my[it+udata->nt*1]-rdata->y[it + udata->nt*1],2.0)*5.0E-1;
  tdata->Jy[1] += rdata->y[it + udata->nt*4]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[2] += rdata->y[it + udata->nt*7]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[3] += rdata->y[it + udata->nt*10]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[4] += rdata->y[it + udata->nt*13]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[5] += rdata->y[it + udata->nt*16]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[6] += rdata->y[it + udata->nt*19]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[7] += rdata->y[it + udata->nt*22]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[8] += rdata->y[it + udata->nt*25]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[9] += rdata->y[it + udata->nt*28]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[10] += rdata->y[it + udata->nt*31]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[11] += rdata->y[it + udata->nt*34]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[12] += rdata->y[it + udata->nt*37]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[13] += rdata->y[it + udata->nt*40]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[14] += rdata->y[it + udata->nt*43]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[15] += rdata->y[it + udata->nt*46]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->Jy[16] += 1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*pow(edata->my[it+udata->nt*1]-rdata->y[it + udata->nt*1],2.0)*-1.0+1.0/tdata->sigmay[1]-rdata->y[it + udata->nt*49]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*5.0E-1;
  tdata->Jy[17] += rdata->y[it + udata->nt*52]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
}
if(!amiIsNaN(edata->my[2* udata->nt+it])){
    iy = 2;
  tdata->Jy[0] += amilog((tdata->sigmay[2]*tdata->sigmay[2])*3.141592653589793*2.0)*5.0E-1+1.0/(tdata->sigmay[2]*tdata->sigmay[2])*pow(edata->my[it+udata->nt*2]-rdata->y[it + udata->nt*2],2.0)*5.0E-1;
  tdata->Jy[1] += rdata->y[it + udata->nt*5]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[2] += rdata->y[it + udata->nt*8]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[3] += rdata->y[it + udata->nt*11]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[4] += rdata->y[it + udata->nt*14]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[5] += rdata->y[it + udata->nt*17]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[6] += rdata->y[it + udata->nt*20]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[7] += rdata->y[it + udata->nt*23]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[8] += rdata->y[it + udata->nt*26]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[9] += rdata->y[it + udata->nt*29]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[10] += rdata->y[it + udata->nt*32]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[11] += rdata->y[it + udata->nt*35]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[12] += rdata->y[it + udata->nt*38]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[13] += rdata->y[it + udata->nt*41]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[14] += rdata->y[it + udata->nt*44]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[15] += rdata->y[it + udata->nt*47]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[16] += rdata->y[it + udata->nt*50]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->Jy[17] += 1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*pow(edata->my[it+udata->nt*2]-rdata->y[it + udata->nt*2],2.0)*-1.0+1.0/tdata->sigmay[2]-rdata->y[it + udata->nt*53]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*5.0E-1;
}
return(status);

}


