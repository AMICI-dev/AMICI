
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_nested_events{

void xdot_model_nested_events(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  xdot[0] = -p[4]*x[0]-p[3]*x[0]*(h[0]-1.0);
}

} // namespace model_model_nested_events

} // namespace amici

