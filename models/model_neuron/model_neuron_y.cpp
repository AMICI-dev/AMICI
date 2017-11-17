
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void y_model_neuron(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  y[0] = x[0];
}

