
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void dxdotdp_model_jakstat_adjoint(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  dxdotdp[0] = -w[0]*x[0];
  dxdotdp[1] = w[0]*x[0];

  } break;

  case 1: {
  dxdotdp[1] = w[1]*-2.0;
  dxdotdp[2] = w[1];

  } break;

  case 2: {
  dxdotdp[2] = -x[2];
  dxdotdp[3] = (k[0]*x[2])/k[1];

  } break;

  case 3: {
  dxdotdp[0] = (k[1]*x[8])/k[0];
  dxdotdp[3] = -x[3];
  dxdotdp[4] = x[3]*2.0-x[4];
  dxdotdp[5] = x[4]-x[5];
  dxdotdp[6] = x[5]-x[6];
  dxdotdp[7] = x[6]-x[7];
  dxdotdp[8] = x[7]-x[8];

  } break;

  case 5: {
  dxdotdp[0] = -p[0]*x[0]*dwdp[0];
  dxdotdp[1] = p[0]*x[0]*dwdp[0];

  } break;

  case 6: {
  dxdotdp[0] = -p[0]*x[0]*dwdp[1];
  dxdotdp[1] = p[0]*x[0]*dwdp[1];

  } break;

  case 7: {
  dxdotdp[0] = -p[0]*x[0]*dwdp[2];
  dxdotdp[1] = p[0]*x[0]*dwdp[2];

  } break;

  case 8: {
  dxdotdp[0] = -p[0]*x[0]*dwdp[3];
  dxdotdp[1] = p[0]*x[0]*dwdp[3];

  } break;

  case 9: {
  dxdotdp[0] = -p[0]*x[0]*dwdp[4];
  dxdotdp[1] = p[0]*x[0]*dwdp[4];

  } break;

}
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

