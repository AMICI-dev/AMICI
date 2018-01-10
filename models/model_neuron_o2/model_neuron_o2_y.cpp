
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void y_model_neuron_o2(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  y[0] = x[0];
  y[1] = x[2];
  y[2] = x[4];
  y[3] = x[6];
  y[4] = x[8];
}

