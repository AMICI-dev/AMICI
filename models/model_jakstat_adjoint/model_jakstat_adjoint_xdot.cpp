
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void xdot_model_jakstat_adjoint(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  xdot[0] = (k[1]*p[3]*x[8]-k[0]*p[0]*w[0]*x[0])/k[0];
  xdot[1] = p[1]*w[1]*-2.0+p[0]*w[0]*x[0];
  xdot[2] = p[1]*w[1]-p[2]*x[2];
  xdot[3] = (k[0]*p[2]*x[2]-k[1]*p[3]*x[3])/k[1];
  xdot[4] = p[3]*(x[3]*2.0-x[4]);
  xdot[5] = p[3]*(x[4]-x[5]);
  xdot[6] = p[3]*(x[5]-x[6]);
  xdot[7] = p[3]*(x[6]-x[7]);
  xdot[8] = p[3]*(x[7]-x[8]);
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

