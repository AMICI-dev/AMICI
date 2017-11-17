
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparse_model_steadystate(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0] = -p[1]*x[1]-p[0]*dwdx[0]*2.0;
  J[1] = -p[1]*x[1]+p[0]*dwdx[0];
  J[2] = p[1]*x[1];
  J[3] = p[2]*2.0-p[1]*x[0];
  J[4] = -p[2]-p[1]*x[0];
  J[5] = p[1]*x[0];
  J[6] = dwdx[1];
  J[7] = dwdx[1];
  J[8] = -k[3]-dwdx[1];
}

