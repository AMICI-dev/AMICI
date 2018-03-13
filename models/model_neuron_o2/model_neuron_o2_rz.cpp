
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void rz_model_neuron_o2(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  rz[0] = x[0]-3.0E1;
  rz[1] = x[2];
  rz[2] = x[4];
  rz[3] = x[6];
  rz[4] = x[8];
}

