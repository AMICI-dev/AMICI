
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void J_model_jakstat_adjoint(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0+0*9] = -p[0]*w[0];
  J[0+8*9] = (k[1]*p[3])/k[0];
  J[1+0*9] = p[0]*w[0];
  J[1+1*9] = p[1]*dwdx[0]*-2.0;
  J[2+1*9] = p[1]*dwdx[0];
  J[2+2*9] = -p[2];
  J[3+2*9] = (k[0]*p[2])/k[1];
  J[3+3*9] = -p[3];
  J[4+3*9] = p[3]*2.0;
  J[4+4*9] = -p[3];
  J[5+4*9] = p[3];
  J[5+5*9] = -p[3];
  J[6+5*9] = p[3];
  J[6+6*9] = -p[3];
  J[7+6*9] = p[3];
  J[7+7*9] = -p[3];
  J[8+7*9] = p[3];
  J[8+8*9] = -p[3];
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

