
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void xBdot_model_neuron_o2(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  xBdot[0] = -xB[0]*(x[0]*(2.0/2.5E1)+5.0)-p[0]*p[1]*xB[1];
  xBdot[1] = xB[0]+p[0]*xB[1];
  xBdot[2] = -p[1]*xB[1]-w[1]*xB[2]-p[0]*p[1]*xB[3]-x[2]*xB[0]*dwdx[1];
  xBdot[3] = xB[1]+xB[2]+p[0]*xB[3];
  xBdot[4] = -p[0]*xB[1]-w[1]*xB[4]-p[0]*p[1]*xB[5]-x[4]*xB[0]*dwdx[1];
  xBdot[5] = xB[4]+p[0]*xB[5];
  xBdot[6] = -w[1]*xB[6]-p[0]*p[1]*xB[7]-x[6]*xB[0]*dwdx[1];
  xBdot[7] = xB[6]+p[0]*xB[7];
  xBdot[8] = -w[1]*xB[8]-p[0]*p[1]*xB[9]-x[8]*xB[0]*dwdx[1];
  xBdot[9] = xB[8]+p[0]*xB[9];
}

