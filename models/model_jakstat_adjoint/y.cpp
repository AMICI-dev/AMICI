
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void y_model_jakstat_adjoint(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  y[0] = p[11]+(p[13]*(x[1]+x[2]*2.0))/p[4];
  y[1] = p[10]+(p[12]*(x[0]+x[1]+x[2]*2.0))/p[4];
  y[2] = amici::spline_pos(t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

