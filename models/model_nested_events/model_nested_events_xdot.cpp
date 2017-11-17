
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void xdot_model_nested_events(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  xdot[0] = -p[4]*x[0]+h[1]*p[3]*x[0];
}

