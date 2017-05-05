
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_steadystate_w.h"

int sJy_model_steadystate(realtype t, int it, realtype *sJy, realtype *s2Jy, realtype *dJydy, realtype *dJydp, realtype *y, realtype *sigma_y, realtype *sy, realtype *dydp, realtype *my, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  sJy[0] -= dJydy[0]*sy[it+nt*0];
  sJy[1] -= dJydy[0]*sy[it+nt*3];
  sJy[2] -= dJydy[0]*sy[it+nt*6];
  sJy[3] -= dJydy[0]*sy[it+nt*9];
  sJy[4] -= dJydy[0]*sy[it+nt*12];
}
if(!amiIsNaN(my[1*nt+it])){
    iy = 1;
  sJy[0] -= dJydy[4]*sy[it+nt*1];
  sJy[1] -= dJydy[4]*sy[it+nt*4];
  sJy[2] -= dJydy[4]*sy[it+nt*7];
  sJy[3] -= dJydy[4]*sy[it+nt*10];
  sJy[4] -= dJydy[4]*sy[it+nt*13];
}
if(!amiIsNaN(my[2*nt+it])){
    iy = 2;
  sJy[0] -= dJydy[8]*sy[it+nt*2];
  sJy[1] -= dJydy[8]*sy[it+nt*5];
  sJy[2] -= dJydy[8]*sy[it+nt*8];
  sJy[3] -= dJydy[8]*sy[it+nt*11];
  sJy[4] -= dJydy[8]*sy[it+nt*14];
}
return(status);

}


