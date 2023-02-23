
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_robertson{

void xdot_model_robertson(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w) {
  xdot[0] = w[0]-dx[0]-p[0]*x[0];
  xdot[1] = -w[0]-dx[1]+p[0]*x[0]-p[2]*(x[1]*x[1]);
  xdot[2] = x[0]+x[1]+x[2]-1.0;
}

} // namespace model_model_robertson

} // namespace amici

