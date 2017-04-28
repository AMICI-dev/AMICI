
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int sJy_model_dirac(realtype t, int it, realtype *sJy, realtype *s2Jy, realtype *dJydy, realtype *dJydp, realtype *y, realtype *sigma_y, realtype *sy, realtype *dydp, realtype *my, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  sJy[0] -= dJydy[0]*sy[it+nt*0];
  sJy[1] -= dJydy[0]*sy[it+nt*1];
  sJy[2] -= dJydy[0]*sy[it+nt*2];
  sJy[3] -= dJydy[0]*sy[it+nt*3];
}
return(status);

}


