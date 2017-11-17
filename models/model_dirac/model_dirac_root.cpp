
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void root_model_dirac(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  root[0] = t-p[1];
  root[1] = -t+p[1];
}

