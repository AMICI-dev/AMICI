
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void dzdx_model_events(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  dzdx[0+1*2] = -1.0/(h[2]+x[2]-p[2]*x[1]+p[1]*x[0]*exp(t*(-1.0/1.0E1))-1.0);
  dzdx[0+2*2] = 1.0/(h[2]+x[2]-p[2]*x[1]+p[1]*x[0]*exp(t*(-1.0/1.0E1))-1.0);
  dzdx[1+0*2] = -1.0/(h[2]+x[2]+p[0]*x[0]*(h[3]-1.0)-1.0);
  dzdx[1+2*2] = 1.0/(h[2]+x[2]+p[0]*x[0]*(h[3]-1.0)-1.0);
}

} // namespace model_model_events

} // namespace amici

