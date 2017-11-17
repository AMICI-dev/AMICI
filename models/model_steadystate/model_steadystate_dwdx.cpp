
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dwdx_model_steadystate(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  dwdx[0] = x[0]*2.0;
  dwdx[1] = p[3];
}

