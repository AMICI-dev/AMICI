
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void qBdot_model_nested_events(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 3: {
  qBdot[0] = -h[1]*x[0]*xB[0];

  } break;

  case 4: {
  qBdot[0] = x[0]*xB[0];

  } break;

}
}

