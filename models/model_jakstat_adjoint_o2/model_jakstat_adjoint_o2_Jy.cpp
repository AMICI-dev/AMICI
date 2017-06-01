
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include "model_jakstat_adjoint_o2_w.h"

int Jy_model_jakstat_adjoint_o2(realtype t, int it, realtype *Jy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0* udata->nt+it])){
    iy = 0;
  Jy[0] += amilog((sigma_y[0]*sigma_y[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[0]*sigma_y[0])*pow(my[it+udata->nt*0]-y[it+udata->nt*0],2.0)*5.0E-1;
  Jy[1] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*3]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[2] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*6]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[3] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*9]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[4] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*12]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[5] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*15]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[6] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*18]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[7] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*21]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[8] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*24]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[9] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*27]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[10] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*30]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[11] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*33]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[12] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*36]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[13] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*39]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[14] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*42]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[15] += 1.0/(sigma_y[0]*sigma_y[0]*sigma_y[0])*pow(my[it+udata->nt*0]-y[it+udata->nt*0],2.0)*-1.0+1.0/sigma_y[0]-1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*45]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*5.0E-1;
  Jy[16] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*48]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
  Jy[17] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+udata->nt*51]*(my[it+udata->nt*0]*2.0-y[it+udata->nt*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[1* udata->nt+it])){
    iy = 1;
  Jy[0] += amilog((sigma_y[1]*sigma_y[1])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[1]*sigma_y[1])*pow(my[it+udata->nt*1]-y[it+udata->nt*1],2.0)*5.0E-1;
  Jy[1] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*4]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[2] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*7]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[3] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*10]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[4] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*13]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[5] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*16]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[6] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*19]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[7] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*22]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[8] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*25]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[9] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*28]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[10] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*31]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[11] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*34]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[12] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*37]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[13] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*40]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[14] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*43]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[15] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*46]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
  Jy[16] += 1.0/(sigma_y[1]*sigma_y[1]*sigma_y[1])*pow(my[it+udata->nt*1]-y[it+udata->nt*1],2.0)*-1.0+1.0/sigma_y[1]-1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*49]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*5.0E-1;
  Jy[17] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+udata->nt*52]*(my[it+udata->nt*1]*2.0-y[it+udata->nt*1]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[2* udata->nt+it])){
    iy = 2;
  Jy[0] += amilog((sigma_y[2]*sigma_y[2])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[2]*sigma_y[2])*pow(my[it+udata->nt*2]-y[it+udata->nt*2],2.0)*5.0E-1;
  Jy[1] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*5]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[2] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*8]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[3] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*11]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[4] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*14]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[5] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*17]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[6] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*20]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[7] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*23]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[8] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*26]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[9] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*29]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[10] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*32]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[11] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*35]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[12] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*38]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[13] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*41]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[14] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*44]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[15] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*47]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[16] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*50]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*-5.0E-1;
  Jy[17] += 1.0/(sigma_y[2]*sigma_y[2]*sigma_y[2])*pow(my[it+udata->nt*2]-y[it+udata->nt*2],2.0)*-1.0+1.0/sigma_y[2]-1.0/(sigma_y[2]*sigma_y[2])*y[it+udata->nt*53]*(my[it+udata->nt*2]*2.0-y[it+udata->nt*2]*2.0)*5.0E-1;
}
return(status);

}


