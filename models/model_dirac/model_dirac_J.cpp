
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_dirac{

void J_model_dirac(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0+0*2] = -p[0];
  J[1+0*2] = p[2];
  J[1+1*2] = -p[3];
}

} // namespace model_model_dirac

} // namespace amici

