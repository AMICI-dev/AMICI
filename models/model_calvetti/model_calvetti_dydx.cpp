
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void dydx_model_calvetti(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  dydx[0+0*6] = 1.0;
  dydx[1+1*6] = 1.0;
  dydx[2+2*6] = 1.0;
  dydx[3+0*6] = (x[3]*2.0)/x[0]+(1.0/(k[0]*k[0])*x[0]*4.0)/k[1]-(1.0/(k[0]*k[0])*x[0]*(x[3]*((k[0]*k[0])*k[1]*1.0/(x[0]*x[0])+(k[2]*k[2])*k[3]*1.0/(x[1]*x[1]))+x[4]*((k[2]*k[2])*k[3]*1.0/(x[1]*x[1])+(k[4]*k[4])*k[5]*1.0/(x[2]*x[2]))+(k[4]*k[4])*k[5]*1.0/(x[2]*x[2])*x[5])*2.0)/k[1];
  dydx[3+1*6] = (1.0/(k[0]*k[0])*(x[0]*x[0])*((k[2]*k[2])*k[3]*1.0/(x[1]*x[1]*x[1])*x[3]*2.0+(k[2]*k[2])*k[3]*1.0/(x[1]*x[1]*x[1])*x[4]*2.0))/k[1];
  dydx[3+2*6] = (1.0/(k[0]*k[0])*(x[0]*x[0])*((k[4]*k[4])*k[5]*1.0/(x[2]*x[2]*x[2])*x[4]*2.0+(k[4]*k[4])*k[5]*1.0/(x[2]*x[2]*x[2])*x[5]*2.0))/k[1];
  dydx[3+3*6] = -(1.0/(k[0]*k[0])*(x[0]*x[0])*((k[0]*k[0])*k[1]*1.0/(x[0]*x[0])+(k[2]*k[2])*k[3]*1.0/(x[1]*x[1])))/k[1];
  dydx[3+4*6] = -(1.0/(k[0]*k[0])*(x[0]*x[0])*((k[2]*k[2])*k[3]*1.0/(x[1]*x[1])+(k[4]*k[4])*k[5]*1.0/(x[2]*x[2])))/k[1];
  dydx[3+5*6] = -(1.0/(k[0]*k[0])*(k[4]*k[4])*k[5]*(x[0]*x[0])*1.0/(x[2]*x[2]))/k[1];
  dydx[4+3*6] = 1.0;
  dydx[5+4*6] = 1.0;
}

} // namespace model_model_calvetti

} // namespace amici

