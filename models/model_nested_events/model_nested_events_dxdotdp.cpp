
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dxdotdp_model_nested_events(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 3: {
  dxdotdp[0] = h[1]*x[0];

  } break;

  case 4: {
  dxdotdp[0] = -x[0];

  } break;

}
}

