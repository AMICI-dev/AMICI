
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void dydx_model_events(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  dydx[0+0*1] = p[3];
  dydx[0+1*1] = p[3];
  dydx[0+2*1] = p[3];
}

} // namespace model_model_events

} // namespace amici

