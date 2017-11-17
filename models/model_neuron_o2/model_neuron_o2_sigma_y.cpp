
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void sigma_y_model_neuron_o2(double *sigmay, const realtype t, const realtype *p, const realtype *k) {
  sigmay[0] = 1.0;
}

