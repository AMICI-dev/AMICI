
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dxdotdp_model_dirac(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  dxdotdp[0] = -x[0];

  } break;

  case 2: {
  dxdotdp[1] = x[0];

  } break;

  case 3: {
  dxdotdp[1] = -x[1];

  } break;

}
}

