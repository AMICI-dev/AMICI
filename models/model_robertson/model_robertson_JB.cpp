
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_robertson{

void JB_model_robertson(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) {
  JB[0+0*3] = cj+p[0];
  JB[0+1*3] = -p[0];
  JB[0+2*3] = -1.0;
  JB[1+0*3] = -dwdx[0];
  JB[1+1*3] = cj+dwdx[0]+p[2]*x[1]*2.0;
  JB[1+2*3] = -1.0;
  JB[2+0*3] = -dwdx[1];
  JB[2+1*3] = dwdx[1];
  JB[2+2*3] = -1.0;
}

} // namespace model_model_robertson

} // namespace amici

