
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_dirac{

void JB_model_dirac(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0+0*2] = p[0];
  JB[0+1*2] = -p[2];
  JB[1+1*2] = p[3];
}

} // namespace model_model_dirac

} // namespace amici

