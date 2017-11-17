
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dydx_model_robertson(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  dydx[0+0*3] = 1.0;
  dydx[1+1*3] = 1.0E4;
  dydx[2+2*3] = 1.0;
}

