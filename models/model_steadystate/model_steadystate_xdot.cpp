
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void xdot_model_steadystate(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  xdot[0] = p[4]+w[0]-p[0]*w[1]*2.0+p[2]*x[1]*2.0-p[1]*x[0]*x[1];
  xdot[1] = w[0]+p[0]*w[1]-p[2]*x[1]-p[1]*x[0]*x[1];
  xdot[2] = -w[0]-k[3]*x[2]+p[1]*x[0]*x[1];
}

