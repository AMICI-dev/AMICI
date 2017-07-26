
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/tdata.h>
#include <include/rdata.h>
#include <include/edata.h>
#include "model_jakstat_adjoint_o2_w.h"

int dJydy_model_jakstat_adjoint_o2(realtype t, int it, N_Vector x, void *user_data, TempData *tdata, const ExpData *edata, ReturnData *rdata) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
memset(tdata->dJydy,0,sizeof(realtype)*udata->ny*udata->nytrue*udata->nJ);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(edata->my[0* udata->nt+it])){
    iy = 0;
  tdata->dJydy[iy+(0+0*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(1+0*18)*udata->nytrue] = rdata->y[it + udata->nt*3]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(1+3*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(2+0*18)*udata->nytrue] = rdata->y[it + udata->nt*6]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(2+6*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(3+0*18)*udata->nytrue] = rdata->y[it + udata->nt*9]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(3+9*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(4+0*18)*udata->nytrue] = rdata->y[it + udata->nt*12]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(4+12*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(5+0*18)*udata->nytrue] = rdata->y[it + udata->nt*15]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(5+15*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(6+0*18)*udata->nytrue] = rdata->y[it + udata->nt*18]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(6+18*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(7+0*18)*udata->nytrue] = rdata->y[it + udata->nt*21]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(7+21*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(8+0*18)*udata->nytrue] = rdata->y[it + udata->nt*24]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(8+24*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(9+0*18)*udata->nytrue] = rdata->y[it + udata->nt*27]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(9+27*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(10+0*18)*udata->nytrue] = rdata->y[it + udata->nt*30]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(10+30*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(11+0*18)*udata->nytrue] = rdata->y[it + udata->nt*33]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(11+33*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(12+0*18)*udata->nytrue] = rdata->y[it + udata->nt*36]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(12+36*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(13+0*18)*udata->nytrue] = rdata->y[it + udata->nt*39]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(13+39*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(14+0*18)*udata->nytrue] = rdata->y[it + udata->nt*42]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(14+42*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(15+0*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*1.0+rdata->y[it + udata->nt*45]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(15+45*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(16+0*18)*udata->nytrue] = rdata->y[it + udata->nt*48]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(16+48*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(17+0*18)*udata->nytrue] = rdata->y[it + udata->nt*51]*1.0/(tdata->sigmay[0]*tdata->sigmay[0])*1.0;
  tdata->dJydy[iy+(17+51*18)*udata->nytrue] = 1.0/(tdata->sigmay[0]*tdata->sigmay[0])*(edata->my[it+udata->nt*0]*2.0-rdata->y[it + udata->nt*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(edata->my[1* udata->nt+it])){
    iy = 1;
  tdata->dJydy[iy+(0+1*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(1+1*18)*udata->nytrue] = rdata->y[it + udata->nt*4]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(1+4*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(2+1*18)*udata->nytrue] = rdata->y[it + udata->nt*7]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(2+7*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(3+1*18)*udata->nytrue] = rdata->y[it + udata->nt*10]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(3+10*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(4+1*18)*udata->nytrue] = rdata->y[it + udata->nt*13]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(4+13*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(5+1*18)*udata->nytrue] = rdata->y[it + udata->nt*16]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(5+16*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(6+1*18)*udata->nytrue] = rdata->y[it + udata->nt*19]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(6+19*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(7+1*18)*udata->nytrue] = rdata->y[it + udata->nt*22]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(7+22*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(8+1*18)*udata->nytrue] = rdata->y[it + udata->nt*25]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(8+25*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(9+1*18)*udata->nytrue] = rdata->y[it + udata->nt*28]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(9+28*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(10+1*18)*udata->nytrue] = rdata->y[it + udata->nt*31]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(10+31*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(11+1*18)*udata->nytrue] = rdata->y[it + udata->nt*34]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(11+34*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(12+1*18)*udata->nytrue] = rdata->y[it + udata->nt*37]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(12+37*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(13+1*18)*udata->nytrue] = rdata->y[it + udata->nt*40]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(13+40*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(14+1*18)*udata->nytrue] = rdata->y[it + udata->nt*43]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(14+43*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(15+1*18)*udata->nytrue] = rdata->y[it + udata->nt*46]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(15+46*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(16+1*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*1.0+rdata->y[it + udata->nt*49]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(16+49*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(17+1*18)*udata->nytrue] = rdata->y[it + udata->nt*52]*1.0/(tdata->sigmay[1]*tdata->sigmay[1])*1.0;
  tdata->dJydy[iy+(17+52*18)*udata->nytrue] = 1.0/(tdata->sigmay[1]*tdata->sigmay[1])*(edata->my[it+udata->nt*1]*2.0-rdata->y[it + udata->nt*1]*2.0)*-5.0E-1;
}
if(!amiIsNaN(edata->my[2* udata->nt+it])){
    iy = 2;
  tdata->dJydy[iy+(0+2*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(1+2*18)*udata->nytrue] = rdata->y[it + udata->nt*5]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(1+5*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(2+2*18)*udata->nytrue] = rdata->y[it + udata->nt*8]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(2+8*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(3+2*18)*udata->nytrue] = rdata->y[it + udata->nt*11]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(3+11*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(4+2*18)*udata->nytrue] = rdata->y[it + udata->nt*14]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(4+14*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(5+2*18)*udata->nytrue] = rdata->y[it + udata->nt*17]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(5+17*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(6+2*18)*udata->nytrue] = rdata->y[it + udata->nt*20]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(6+20*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(7+2*18)*udata->nytrue] = rdata->y[it + udata->nt*23]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(7+23*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(8+2*18)*udata->nytrue] = rdata->y[it + udata->nt*26]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(8+26*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(9+2*18)*udata->nytrue] = rdata->y[it + udata->nt*29]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(9+29*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(10+2*18)*udata->nytrue] = rdata->y[it + udata->nt*32]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(10+32*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(11+2*18)*udata->nytrue] = rdata->y[it + udata->nt*35]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(11+35*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(12+2*18)*udata->nytrue] = rdata->y[it + udata->nt*38]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(12+38*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(13+2*18)*udata->nytrue] = rdata->y[it + udata->nt*41]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(13+41*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(14+2*18)*udata->nytrue] = rdata->y[it + udata->nt*44]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(14+44*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(15+2*18)*udata->nytrue] = rdata->y[it + udata->nt*47]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(15+47*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(16+2*18)*udata->nytrue] = rdata->y[it + udata->nt*50]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(16+50*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
  tdata->dJydy[iy+(17+2*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*1.0+rdata->y[it + udata->nt*53]*1.0/(tdata->sigmay[2]*tdata->sigmay[2])*1.0;
  tdata->dJydy[iy+(17+53*18)*udata->nytrue] = 1.0/(tdata->sigmay[2]*tdata->sigmay[2])*(edata->my[it+udata->nt*2]*2.0-rdata->y[it + udata->nt*2]*2.0)*-5.0E-1;
}
return(status);

}


