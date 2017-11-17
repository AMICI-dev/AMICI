
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void qBdot_model_steadystate(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  qBdot[0] = w[1]*xB[0]*2.0-w[1]*xB[1];

  } break;

  case 1: {
  qBdot[0] = x[0]*x[1]*xB[0]+x[0]*x[1]*xB[1]-x[0]*x[1]*xB[2];

  } break;

  case 2: {
  qBdot[0] = x[1]*xB[0]*-2.0+x[1]*xB[1];

  } break;

  case 3: {
  qBdot[0] = -xB[0]*dwdp[0]-xB[1]*dwdp[0]+xB[2]*dwdp[0];

  } break;

  case 4: {
  qBdot[0] = -xB[0];

  } break;

}
}

