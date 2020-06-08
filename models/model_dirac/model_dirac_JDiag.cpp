
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_dirac{

void JDiag_model_dirac(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JDiag[0+0*2] = -p[0];
  JDiag[1+0*2] = -p[3];
}

} // namespace model_model_dirac

} // namespace amici

