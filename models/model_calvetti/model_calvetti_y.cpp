
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void y_model_calvetti(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];
  y[3] = (1.0/(k[0]*k[0])*(x[0]*x[0])*2.0)/k[1]-(1.0/(k[0]*k[0])*(x[0]*x[0])*(x[3]*((k[0]*k[0])*k[1]*1.0/(x[0]*x[0])+(k[2]*k[2])*k[3]*1.0/(x[1]*x[1]))+x[4]*((k[2]*k[2])*k[3]*1.0/(x[1]*x[1])+(k[4]*k[4])*k[5]*1.0/(x[2]*x[2]))+(k[4]*k[4])*k[5]*1.0/(x[2]*x[2])*x[5]))/k[1];
  y[4] = x[3];
  y[5] = x[4];
}

} // namespace model_model_calvetti

} // namespace amici

