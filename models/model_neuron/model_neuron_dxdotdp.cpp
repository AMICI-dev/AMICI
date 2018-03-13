
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dxdotdp_model_neuron(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  dxdotdp[1] = -x[1]+p[1]*x[0];

  } break;

  case 1: {
  dxdotdp[1] = p[0]*x[0];

  } break;

}
}

