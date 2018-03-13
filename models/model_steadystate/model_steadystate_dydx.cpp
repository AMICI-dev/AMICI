
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dydx_model_steadystate(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  dydx[0+0*3] = 1.0;
  dydx[1+1*3] = 1.0;
  dydx[2+2*3] = 1.0;
}

