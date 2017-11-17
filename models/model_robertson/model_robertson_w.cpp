
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void w_model_robertson(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  w[0] = p[1]*x[1]*x[2];
}

