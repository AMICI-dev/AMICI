
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void dxdotdp_model_neuron_o2(realtype *dxdotdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  dxdotdp[1] = -x[1]+p[1]*x[0];
  dxdotdp[3] = -x[3]+p[1]*x[2];
  dxdotdp[5] = x[0]-x[5]+p[1]*x[4];
  dxdotdp[7] = -x[7]+p[1]*x[6];
  dxdotdp[9] = -x[9]+p[1]*x[8];

  } break;

  case 1: {
  dxdotdp[1] = p[0]*x[0];
  dxdotdp[3] = x[0]+p[0]*x[2];
  dxdotdp[5] = p[0]*x[4];
  dxdotdp[7] = p[0]*x[6];
  dxdotdp[9] = p[0]*x[8];

  } break;

}
}

} // namespace model_model_neuron_o2

} // namespace amici

