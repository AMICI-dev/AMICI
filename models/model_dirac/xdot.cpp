
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_dirac{

void xdot_model_dirac(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  xdot[0] = -p[0]*x[0];
  xdot[1] = p[2]*x[0]-p[3]*x[1];
}

} // namespace model_model_dirac

} // namespace amici

