
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparseB_model_jakstat_adjoint(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0] = p[0]*w[0];
  JB[1] = -(k[1]*p[3])/k[0];
  JB[2] = -p[0]*w[0];
  JB[3] = p[1]*dwdx[0]*2.0;
  JB[4] = -p[1]*dwdx[0];
  JB[5] = p[2];
  JB[6] = -(k[0]*p[2])/k[1];
  JB[7] = p[3];
  JB[8] = p[3]*-2.0;
  JB[9] = p[3];
  JB[10] = -p[3];
  JB[11] = p[3];
  JB[12] = -p[3];
  JB[13] = p[3];
  JB[14] = -p[3];
  JB[15] = p[3];
  JB[16] = -p[3];
  JB[17] = p[3];
}

