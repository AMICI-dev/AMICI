
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void qBdot_model_dirac(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  qBdot[0] = x[0]*xB[0];

  } break;

  case 2: {
  qBdot[0] = -x[0]*xB[1];

  } break;

  case 3: {
  qBdot[0] = x[1]*xB[1];

  } break;

}
}

