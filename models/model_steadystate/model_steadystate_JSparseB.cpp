
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparseB_model_steadystate(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0] = p[1]*x[1]+p[0]*dwdx[0]*2.0;
  JB[1] = p[2]*-2.0+p[1]*x[0];
  JB[2] = -dwdx[1];
  JB[3] = p[1]*x[1]-p[0]*dwdx[0];
  JB[4] = p[2]+p[1]*x[0];
  JB[5] = -dwdx[1];
  JB[6] = -p[1]*x[1];
  JB[7] = -p[1]*x[0];
  JB[8] = k[3]+dwdx[1];
}

