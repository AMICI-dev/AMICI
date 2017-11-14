
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JB_model_events(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0+0*3] = h[3]*p[0];
  JB[0+1*3] = -p[1]*exp(t*(-1.0/1.0E1));
  JB[1+1*3] = p[2];
  JB[2+2*3] = 1.0;
}

