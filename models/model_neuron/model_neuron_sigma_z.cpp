
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void sigma_z_model_neuron(double *sigmaz, const realtype t, const realtype *p, const realtype *k) {
  sigmaz[0] = 1.0;
}

