
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void xBdot_model_events(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  xBdot[0] = -p[1]*xB[1]*exp(t*(-1.0/1.0E1))+h[3]*p[0]*xB[0];
  xBdot[1] = p[2]*xB[1];
  xBdot[2] = xB[2];
}

