
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dydx_model_neuron(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  dydx[0+0*1] = 1.0;
}

