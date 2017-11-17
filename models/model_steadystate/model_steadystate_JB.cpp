
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JB_model_steadystate(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0+0*3] = p[1]*x[1]+p[0]*dwdx[0]*2.0;
  JB[0+1*3] = p[1]*x[1]-p[0]*dwdx[0];
  JB[0+2*3] = -p[1]*x[1];
  JB[1+0*3] = p[2]*-2.0+p[1]*x[0];
  JB[1+1*3] = p[2]+p[1]*x[0];
  JB[1+2*3] = -p[1]*x[0];
  JB[2+0*3] = -dwdx[1];
  JB[2+1*3] = -dwdx[1];
  JB[2+2*3] = k[3]+dwdx[1];
}

