
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_robertson{

void JDiag_model_robertson(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) {
  JDiag[0+0*3] = -cj-p[0];
  JDiag[1+0*3] = -cj-dwdx[0]-p[2]*x[1]*2.0;
  JDiag[2+0*3] = 1.0;
}

} // namespace model_model_robertson

} // namespace amici

