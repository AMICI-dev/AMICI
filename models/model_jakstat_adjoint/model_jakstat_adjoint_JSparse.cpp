
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparse_model_jakstat_adjoint(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0] = -p[0]*w[0];
  J[1] = p[0]*w[0];
  J[2] = p[1]*dwdx[0]*-2.0;
  J[3] = p[1]*dwdx[0];
  J[4] = -p[2];
  J[5] = (k[0]*p[2])/k[1];
  J[6] = -p[3];
  J[7] = p[3]*2.0;
  J[8] = -p[3];
  J[9] = p[3];
  J[10] = -p[3];
  J[11] = p[3];
  J[12] = -p[3];
  J[13] = p[3];
  J[14] = -p[3];
  J[15] = p[3];
  J[16] = (k[1]*p[3])/k[0];
  J[17] = -p[3];
}

