
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dwdx_model_jakstat_adjoint_o2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  dwdx[0] = x[1]*2.0;
  dwdx[1] = 2.0;
}

