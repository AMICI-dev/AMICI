
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dydx_model_neuron_o2(double *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  dydx[0+0*5] = 1.0;
  dydx[1+2*5] = 1.0;
  dydx[2+4*5] = 1.0;
  dydx[3+6*5] = 1.0;
  dydx[4+8*5] = 1.0;
}

