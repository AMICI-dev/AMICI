
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void xBdot_model_robertson(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) {
  xBdot[0] = -xB[2]+dxB[0]+p[0]*xB[0]-p[0]*xB[1];
  xBdot[1] = -xB[2]+dxB[1]-xB[0]*dwdx[0]+xB[1]*(dwdx[0]+p[2]*x[1]*2.0);
  xBdot[2] = -xB[2]-xB[0]*dwdx[1]+xB[1]*dwdx[1];
}

