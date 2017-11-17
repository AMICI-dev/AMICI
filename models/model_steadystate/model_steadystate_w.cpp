
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void w_model_steadystate(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  w[0] = p[3]*x[2];
  w[1] = x[0]*x[0];
}

