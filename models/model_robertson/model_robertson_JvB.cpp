
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void JvB_model_robertson(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *vB, const realtype *w, const realtype *dwdx) {
  JvB[0] = -vB[2]+vB[0]*(cj+p[0])-p[0]*vB[1];
  JvB[1] = -vB[2]-vB[0]*dwdx[0]+vB[1]*(cj+dwdx[0]+p[2]*x[1]*2.0);
  JvB[2] = -vB[2]-vB[0]*dwdx[1]+vB[1]*dwdx[1];
}

