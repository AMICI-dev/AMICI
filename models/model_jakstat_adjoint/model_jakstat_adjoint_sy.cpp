
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_w.h"

int sy_model_jakstat_adjoint(realtype t, int it, realtype *sy, realtype *dydx, realtype *dydp, N_Vector *sx, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *sx_tmp;
int ip;
for(ip = 0; ip<nplist; ip++) {
sx_tmp = N_VGetArrayPointer(sx[ip]);
switch (plist[ip]) {
  case 0: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 1: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 2: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 3: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 4: {
  sy[it+nt*((0)+ip*3)] = dydp[12]+dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydp[13]+dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 5: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];
  sy[it+nt*((2)+ip*3)] = dydp[17];

  } break;

  case 6: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];
  sy[it+nt*((2)+ip*3)] = dydp[20];

  } break;

  case 7: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];
  sy[it+nt*((2)+ip*3)] = dydp[23];

  } break;

  case 8: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];
  sy[it+nt*((2)+ip*3)] = dydp[26];

  } break;

  case 9: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];
  sy[it+nt*((2)+ip*3)] = dydp[29];

  } break;

  case 10: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydp[31]+dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 11: {
  sy[it+nt*((0)+ip*3)] = dydp[33]+dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 12: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydp[37]+dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 13: {
  sy[it+nt*((0)+ip*3)] = dydp[39]+dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 14: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 15: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

  case 16: {
  sy[it+nt*((0)+ip*3)] = dydx[3]*sx_tmp[1]+dydx[6]*sx_tmp[2];
  sy[it+nt*((1)+ip*3)] = dydx[1]*sx_tmp[0]+dydx[4]*sx_tmp[1]+dydx[7]*sx_tmp[2];

  } break;

}
}
return(status);

}


