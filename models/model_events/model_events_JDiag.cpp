
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void JDiag_model_events(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JDiag[0+0*3] = p[0]*(h[3]-1.0);
  JDiag[1+0*3] = -p[2];
  JDiag[2+0*3] = -1.0;
}

} // namespace model_model_events

} // namespace amici

