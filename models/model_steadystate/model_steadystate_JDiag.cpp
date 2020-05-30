
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_steadystate{

void JDiag_model_steadystate(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JDiag[0+0*3] = -p[1]*x[1]-p[0]*dwdx[0]*2.0;
  JDiag[1+0*3] = -p[2]-p[1]*x[0];
  JDiag[2+0*3] = -k[3]-dwdx[1];
}

} // namespace model_model_steadystate

} // namespace amici 

