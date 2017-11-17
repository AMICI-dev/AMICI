
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparseB_model_robertson(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) {
  JB[0] = cj+p[0];
  JB[1] = -dwdx[0];
  JB[2] = -dwdx[1];
  JB[3] = -p[0];
  JB[4] = cj+dwdx[0]+p[2]*x[1]*2.0;
  JB[5] = dwdx[1];
  JB[6] = -1.0;
  JB[7] = -1.0;
  JB[8] = -1.0;
}

