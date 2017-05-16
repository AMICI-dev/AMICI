
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int sJy_model_jakstat_adjoint(realtype t, int it, realtype *sJy, realtype *s2Jy, realtype *dJydy, realtype *dJydp, realtype *y, realtype *sigma_y, realtype *sy, realtype *dydp, realtype *my, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  sJy[0] -= dJydy[0]*sy[it+nt*0];
  sJy[1] -= dJydy[0]*sy[it+nt*3];
  sJy[2] -= dJydy[0]*sy[it+nt*6];
  sJy[3] -= dJydy[0]*sy[it+nt*9];
  sJy[4] -= dJydp[12]-dJydy[0]*(dydp[12]-sy[it+nt*12]);
  sJy[5] -= dJydy[0]*sy[it+nt*15];
  sJy[6] -= dJydy[0]*sy[it+nt*18];
  sJy[7] -= dJydy[0]*sy[it+nt*21];
  sJy[8] -= dJydy[0]*sy[it+nt*24];
  sJy[9] -= dJydy[0]*sy[it+nt*27];
  sJy[10] -= dJydy[0]*sy[it+nt*30];
  sJy[11] -= dJydp[33]-dJydy[0]*(dydp[33]-sy[it+nt*33]);
  sJy[12] -= dJydy[0]*sy[it+nt*36];
  sJy[13] -= dJydp[39]-dJydy[0]*(dydp[39]-sy[it+nt*39]);
  sJy[14] -= dJydp[42]+dJydy[0]*sy[it+nt*42];
  sJy[15] -= dJydy[0]*sy[it+nt*45];
  sJy[16] -= dJydy[0]*sy[it+nt*48];
}
if(!amiIsNaN(my[1*nt+it])){
    iy = 1;
  sJy[0] -= dJydy[4]*sy[it+nt*1];
  sJy[1] -= dJydy[4]*sy[it+nt*4];
  sJy[2] -= dJydy[4]*sy[it+nt*7];
  sJy[3] -= dJydy[4]*sy[it+nt*10];
  sJy[4] -= dJydp[13]-dJydy[4]*(dydp[13]-sy[it+nt*13]);
  sJy[5] -= dJydy[4]*sy[it+nt*16];
  sJy[6] -= dJydy[4]*sy[it+nt*19];
  sJy[7] -= dJydy[4]*sy[it+nt*22];
  sJy[8] -= dJydy[4]*sy[it+nt*25];
  sJy[9] -= dJydy[4]*sy[it+nt*28];
  sJy[10] -= dJydp[31]-dJydy[4]*(dydp[31]-sy[it+nt*31]);
  sJy[11] -= dJydy[4]*sy[it+nt*34];
  sJy[12] -= dJydp[37]-dJydy[4]*(dydp[37]-sy[it+nt*37]);
  sJy[13] -= dJydy[4]*sy[it+nt*40];
  sJy[14] -= dJydy[4]*sy[it+nt*43];
  sJy[15] -= dJydp[46]+dJydy[4]*sy[it+nt*46];
  sJy[16] -= dJydy[4]*sy[it+nt*49];
}
if(!amiIsNaN(my[2*nt+it])){
    iy = 2;
  sJy[5] -= dJydp[17]-dJydy[8]*(dydp[17]-sy[it+nt*17]);
  sJy[6] -= dJydp[20]-dJydy[8]*(dydp[20]-sy[it+nt*20]);
  sJy[7] -= dJydp[23]-dJydy[8]*(dydp[23]-sy[it+nt*23]);
  sJy[8] -= dJydp[26]-dJydy[8]*(dydp[26]-sy[it+nt*26]);
  sJy[9] -= dJydp[29]-dJydy[8]*(dydp[29]-sy[it+nt*29]);
  sJy[16] -= dJydp[50];
}
return(status);

}


