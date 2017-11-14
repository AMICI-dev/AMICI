
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JSparse_model_events(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0] = -h[3]*p[0];
  J[1] = p[1]*exp(t*(-1.0/1.0E1));
  J[2] = -p[2];
  J[3] = -1.0;
}

