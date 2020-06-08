
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_steadystate{

void dxdotdp_model_steadystate(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  dxdotdp[0] = w[1]*-2.0;
  dxdotdp[1] = w[1];

  } break;

  case 1: {
  dxdotdp[0] = -x[0]*x[1];
  dxdotdp[1] = -x[0]*x[1];
  dxdotdp[2] = x[0]*x[1];

  } break;

  case 2: {
  dxdotdp[0] = x[1]*2.0;
  dxdotdp[1] = -x[1];

  } break;

  case 3: {
  dxdotdp[0] = dwdp[0];
  dxdotdp[1] = dwdp[0];
  dxdotdp[2] = -dwdp[0];

  } break;

  case 4: {
  dxdotdp[0] = 1.0;

  } break;

}
}

} // namespace model_model_steadystate

} // namespace amici

