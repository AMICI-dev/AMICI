
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparseB_model_events(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0] = h[3]*p[0];
  JB[1] = -p[1]*exp(t*(-1.0/1.0E1));
  JB[2] = p[2];
  JB[3] = 1.0;
}

