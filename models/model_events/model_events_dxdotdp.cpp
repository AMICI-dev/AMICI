
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void dxdotdp_model_events(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  dxdotdp[0] = x[0]*(h[3]-1.0);

  } break;

  case 1: {
  dxdotdp[1] = x[0]*exp(t*(-1.0/1.0E1));

  } break;

  case 2: {
  dxdotdp[1] = -x[1];

  } break;

}
}

} // namespace model_model_events

} // namespace amici

