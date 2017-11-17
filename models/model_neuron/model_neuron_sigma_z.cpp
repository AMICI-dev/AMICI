
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void sigma_z_model_neuron(double *sigmaz, const realtype t, const realtype *p, const realtype *k) {
  sigmaz[0] = 1.0;
}

