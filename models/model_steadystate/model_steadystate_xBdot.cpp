
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void xBdot_model_steadystate(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  xBdot[0] = xB[0]*(p[1]*x[1]+p[0]*dwdx[0]*2.0)+xB[1]*(p[1]*x[1]-p[0]*dwdx[0])-p[1]*x[1]*xB[2];
  xBdot[1] = -xB[0]*(p[2]*2.0-p[1]*x[0])+xB[1]*(p[2]+p[1]*x[0])-p[1]*x[0]*xB[2];
  xBdot[2] = xB[2]*(k[3]+dwdx[1])-xB[0]*dwdx[1]-xB[1]*dwdx[1];
}

