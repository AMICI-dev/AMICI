
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparse_model_robertson(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) {
  J[0] = -cj-p[0];
  J[1] = p[0];
  J[2] = 1.0;
  J[3] = dwdx[0];
  J[4] = -cj-dwdx[0]-p[2]*x[1]*2.0;
  J[5] = 1.0;
  J[6] = dwdx[1];
  J[7] = -dwdx[1];
  J[8] = 1.0;
}

