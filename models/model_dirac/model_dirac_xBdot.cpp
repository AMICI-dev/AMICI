
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void xBdot_model_dirac(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  xBdot[0] = p[0]*xB[0]-p[2]*xB[1];
  xBdot[1] = p[3]*xB[1];
}

