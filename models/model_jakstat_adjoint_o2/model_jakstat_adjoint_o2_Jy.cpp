
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int Jy_model_jakstat_adjoint_o2(realtype t, int it, realtype *Jy, realtype *y, N_Vector x, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  Jy[0] += amilog((sigma_y[0]*sigma_y[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[0]*sigma_y[0])*pow(my[it+nt*0]-y[it+nt*0],2.0)*5.0E-1;
  Jy[1] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*3]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[2] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*6]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[3] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*9]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[4] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*12]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[5] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*15]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[6] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*18]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[7] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*21]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[8] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*24]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[9] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*27]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[10] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*30]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[11] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*33]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[12] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*36]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[13] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*39]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[14] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*42]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[15] += 1.0/(sigma_y[0]*sigma_y[0]*sigma_y[0])*pow(my[it+nt*0]-y[it+nt*0],2.0)*-1.0+1.0/sigma_y[0]-1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*45]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*5.0E-1;
  Jy[16] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*48]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  Jy[17] += 1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*51]*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[1*nt+it])){
    iy = 1;
  Jy[0] += amilog((sigma_y[1]*sigma_y[1])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[1]*sigma_y[1])*pow(my[it+nt*1]-y[it+nt*1],2.0)*5.0E-1;
  Jy[1] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*4]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[2] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*7]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[3] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*10]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[4] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*13]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[5] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*16]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[6] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*19]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[7] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*22]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[8] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*25]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[9] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*28]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[10] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*31]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[11] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*34]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[12] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*37]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[13] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*40]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[14] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*43]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[15] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*46]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  Jy[16] += 1.0/(sigma_y[1]*sigma_y[1]*sigma_y[1])*pow(my[it+nt*1]-y[it+nt*1],2.0)*-1.0+1.0/sigma_y[1]-1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*49]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*5.0E-1;
  Jy[17] += 1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*52]*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
}
if(!amiIsNaN(my[2*nt+it])){
    iy = 2;
  Jy[0] += amilog((sigma_y[2]*sigma_y[2])*3.141592653589793*2.0)*5.0E-1+1.0/(sigma_y[2]*sigma_y[2])*pow(my[it+nt*2]-y[it+nt*2],2.0)*5.0E-1;
  Jy[1] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*5]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[2] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*8]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[3] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*11]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[4] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*14]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[5] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*17]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[6] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*20]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[7] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*23]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[8] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*26]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[9] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*29]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[10] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*32]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[11] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*35]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[12] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*38]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[13] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*41]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[14] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*44]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[15] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*47]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[16] += 1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*50]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*-5.0E-1;
  Jy[17] += 1.0/(sigma_y[2]*sigma_y[2]*sigma_y[2])*pow(my[it+nt*2]-y[it+nt*2],2.0)*-1.0+1.0/sigma_y[2]-1.0/(sigma_y[2]*sigma_y[2])*y[it+nt*53]*(my[it+nt*2]*2.0-y[it+nt*2]*2.0)*5.0E-1;
}
return(status);

}


