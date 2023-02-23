
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_steadystate{

void y_model_steadystate(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  y[0] = x[0];
  y[1] = x[1];
  y[2] = x[2];
}

} // namespace model_model_steadystate

} // namespace amici

