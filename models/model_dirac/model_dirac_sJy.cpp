
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_dirac_w.h"

int sJy_model_dirac(realtype t, int it, realtype *sJy, realtype *dJydy, realtype *dJydp, realtype *sy, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
  sJy[0+0*1] = dJydy[0]*sy[0];
  sJy[0+1*1] = dJydy[0]*sy[1];
  sJy[0+2*1] = dJydy[0]*sy[2];
  sJy[0+3*1] = dJydy[0]*sy[3];
return(status);

}


