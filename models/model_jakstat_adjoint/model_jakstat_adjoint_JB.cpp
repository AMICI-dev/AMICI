
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void JB_model_jakstat_adjoint(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0+0*9] = p[0]*w[0];
  JB[0+1*9] = -p[0]*w[0];
  JB[1+1*9] = p[1]*dwdx[0]*2.0;
  JB[1+2*9] = -p[1]*dwdx[0];
  JB[2+2*9] = p[2];
  JB[2+3*9] = -(k[0]*p[2])/k[1];
  JB[3+3*9] = p[3];
  JB[3+4*9] = p[3]*-2.0;
  JB[4+4*9] = p[3];
  JB[4+5*9] = -p[3];
  JB[5+5*9] = p[3];
  JB[5+6*9] = -p[3];
  JB[6+6*9] = p[3];
  JB[6+7*9] = -p[3];
  JB[7+7*9] = p[3];
  JB[7+8*9] = -p[3];
  JB[8+0*9] = -(k[1]*p[3])/k[0];
  JB[8+8*9] = p[3];
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

