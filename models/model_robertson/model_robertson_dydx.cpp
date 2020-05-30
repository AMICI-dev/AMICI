
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_robertson{

void dydx_model_robertson(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  dydx[0+0*3] = 1.0;
  dydx[1+1*3] = 1.0E4;
  dydx[2+2*3] = 1.0;
}

} // namespace model_model_robertson

} // namespace amici

