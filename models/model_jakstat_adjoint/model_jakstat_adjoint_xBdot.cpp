
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void xBdot_model_jakstat_adjoint(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  xBdot[0] = p[0]*w[0]*xB[0]-p[0]*w[0]*xB[1];
  xBdot[1] = p[1]*xB[1]*dwdx[0]*2.0-p[1]*xB[2]*dwdx[0];
  xBdot[2] = p[2]*xB[2]-(k[0]*p[2]*xB[3])/k[1];
  xBdot[3] = p[3]*xB[3]-p[3]*xB[4]*2.0;
  xBdot[4] = p[3]*xB[4]-p[3]*xB[5];
  xBdot[5] = p[3]*xB[5]-p[3]*xB[6];
  xBdot[6] = p[3]*xB[6]-p[3]*xB[7];
  xBdot[7] = p[3]*xB[7]-p[3]*xB[8];
  xBdot[8] = p[3]*xB[8]-(k[1]*p[3]*xB[0])/k[0];
}

