
#include <include/symbolic_functions.h>
#include <include/amici.h>
#include <include/amici_model.h>
#include <string.h>
#include <include/tdata.h>
#include <include/udata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_jakstat_adjoint_o2_w.h"

using namespace amici;

void dJydsigma_model_jakstat_adjoint_o2(realtype t, int it, N_Vector x, amici::TempData *tdata, const amici::ExpData *edata, amici::ReturnData *rdata) {
Model *model = (Model*) tdata->model;
UserData *udata = (UserData*) tdata->udata;
realtype *x_tmp = nullptr;
if(x)
    x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJydsigma,0,sizeof(realtype)*model->nytrue*model->ny*model->nJ);
w_model_jakstat_adjoint_o2(t,x,NULL,tdata);
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->dJydsigma[iy+(0+0*18)*model->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*pow(edata->my[it+udata->nt*0]-rdata->y[it + udata->nt*0],2.0)*-1.0+1.0/tdata->sigmay[0];
  tdata->dJydsigma[iy+(1+0*18)*model->nytrue] = rdata->y[it + udata->nt*3]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(2+0*18)*model->nytrue] = rdata->y[it + udata->nt*6]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(3+0*18)*model->nytrue] = rdata->y[it + udata->nt*9]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(4+0*18)*model->nytrue] = rdata->y[it + udata->nt*12]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(5+0*18)*model->nytrue] = rdata->y[it + udata->nt*15]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(6+0*18)*model->nytrue] = rdata->y[it + udata->nt*18]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(7+0*18)*model->nytrue] = rdata->y[it + udata->nt*21]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(8+0*18)*model->nytrue] = rdata->y[it + udata->nt*24]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(9+0*18)*model->nytrue] = rdata->y[it + udata->nt*27]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(10+0*18)*model->nytrue] = rdata->y[it + udata->nt*30]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(11+0*18)*model->nytrue] = rdata->y[it + udata->nt*33]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(12+0*18)*model->nytrue] = rdata->y[it + udata->nt*36]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(13+0*18)*model->nytrue] = rdata->y[it + udata->nt*39]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(14+0*18)*model->nytrue] = rdata->y[it + udata->nt*42]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(15+0*18)*model->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*pow(edata->my[it+udata->nt*0]-rdata->y[it + udata->nt*0],2.0)*3.0-1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0+rdata->y[it + udata->nt*45]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(16+0*18)*model->nytrue] = rdata->y[it + udata->nt*48]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
  tdata->dJydsigma[iy+(17+0*18)*model->nytrue] = rdata->y[it + udata->nt*51]*1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0;
}
if(!amiIsNaN(edata->my[1* udata->nt+it])){
    iy = 1;
  tdata->dJydsigma[iy+(0+1*18)*model->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*pow(edata->my[it+udata->nt*1]-rdata->y[it + udata->nt*1],2.0)*-1.0+1.0/tdata->sigmay[1];
  tdata->dJydsigma[iy+(1+1*18)*model->nytrue] = rdata->y[it + udata->nt*4]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(2+1*18)*model->nytrue] = rdata->y[it + udata->nt*7]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(3+1*18)*model->nytrue] = rdata->y[it + udata->nt*10]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(4+1*18)*model->nytrue] = rdata->y[it + udata->nt*13]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(5+1*18)*model->nytrue] = rdata->y[it + udata->nt*16]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(6+1*18)*model->nytrue] = rdata->y[it + udata->nt*19]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(7+1*18)*model->nytrue] = rdata->y[it + udata->nt*22]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(8+1*18)*model->nytrue] = rdata->y[it + udata->nt*25]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(9+1*18)*model->nytrue] = rdata->y[it + udata->nt*28]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(10+1*18)*model->nytrue] = rdata->y[it + udata->nt*31]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(11+1*18)*model->nytrue] = rdata->y[it + udata->nt*34]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(12+1*18)*model->nytrue] = rdata->y[it + udata->nt*37]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(13+1*18)*model->nytrue] = rdata->y[it + udata->nt*40]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(14+1*18)*model->nytrue] = rdata->y[it + udata->nt*43]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(15+1*18)*model->nytrue] = rdata->y[it + udata->nt*46]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(16+1*18)*model->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*pow(edata->my[it+udata->nt*1]-rdata->y[it + udata->nt*1],2.0)*3.0-1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0+rdata->y[it + udata->nt*49]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
  tdata->dJydsigma[iy+(17+1*18)*model->nytrue] = rdata->y[it + udata->nt*52]*1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0;
}
if(!amiIsNaN(edata->my[2* udata->nt+it])){
    iy = 2;
  tdata->dJydsigma[iy+(0+2*18)*model->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*pow(edata->my[it+udata->nt*2]-rdata->y[it + udata->nt*2],2.0)*-1.0+1.0/tdata->sigmay[2];
  tdata->dJydsigma[iy+(1+2*18)*model->nytrue] = rdata->y[it + udata->nt*5]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(2+2*18)*model->nytrue] = rdata->y[it + udata->nt*8]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(3+2*18)*model->nytrue] = rdata->y[it + udata->nt*11]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(4+2*18)*model->nytrue] = rdata->y[it + udata->nt*14]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(5+2*18)*model->nytrue] = rdata->y[it + udata->nt*17]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(6+2*18)*model->nytrue] = rdata->y[it + udata->nt*20]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(7+2*18)*model->nytrue] = rdata->y[it + udata->nt*23]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(8+2*18)*model->nytrue] = rdata->y[it + udata->nt*26]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(9+2*18)*model->nytrue] = rdata->y[it + udata->nt*29]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(10+2*18)*model->nytrue] = rdata->y[it + udata->nt*32]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(11+2*18)*model->nytrue] = rdata->y[it + udata->nt*35]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(12+2*18)*model->nytrue] = rdata->y[it + udata->nt*38]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(13+2*18)*model->nytrue] = rdata->y[it + udata->nt*41]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(14+2*18)*model->nytrue] = rdata->y[it + udata->nt*44]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(15+2*18)*model->nytrue] = rdata->y[it + udata->nt*47]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(16+2*18)*model->nytrue] = rdata->y[it + udata->nt*50]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
  tdata->dJydsigma[iy+(17+2*18)*model->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*pow(edata->my[it+udata->nt*2]-rdata->y[it + udata->nt*2],2.0)*3.0-1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0+rdata->y[it + udata->nt*53]*1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0;
}
return;

}


