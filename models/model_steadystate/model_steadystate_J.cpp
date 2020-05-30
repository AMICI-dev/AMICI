
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_steadystate{

void J_model_steadystate(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0+0*3] = -p[1]*x[1]-p[0]*dwdx[0]*2.0;
  J[0+1*3] = p[2]*2.0-p[1]*x[0];
  J[0+2*3] = dwdx[1];
  J[1+0*3] = -p[1]*x[1]+p[0]*dwdx[0];
  J[1+1*3] = -p[2]-p[1]*x[0];
  J[1+2*3] = dwdx[1];
  J[2+0*3] = p[1]*x[1];
  J[2+1*3] = p[1]*x[0];
  J[2+2*3] = -k[3]-dwdx[1];
}

} // namespace model_model_steadystate

} // namespace amici

