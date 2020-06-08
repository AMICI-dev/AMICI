
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void xdot_model_neuron_o2(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  xdot[0] = k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2;
  xdot[1] = -p[0]*(x[1]-p[1]*x[0]);
  xdot[2] = -x[3]+w[1]*x[2];
  xdot[3] = -x[1]+p[1]*x[0]-p[0]*x[3]+p[0]*p[1]*x[2];
  xdot[4] = -x[5]+w[1]*x[4];
  xdot[5] = p[0]*x[0]-p[0]*x[5]+p[0]*p[1]*x[4];
  xdot[6] = -x[7]+w[1]*x[6];
  xdot[7] = -p[0]*x[7]+p[0]*p[1]*x[6];
  xdot[8] = -x[9]+w[1]*x[8];
  xdot[9] = -p[0]*x[9]+p[0]*p[1]*x[8];
}

} // namespace model_model_neuron_o2

} // namespace amici

