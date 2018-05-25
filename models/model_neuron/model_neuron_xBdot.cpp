
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void xBdot_model_neuron(realtype *xBdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  xBdot[0] = -xB[0]*(x[0]*(2.0/2.5E1)+5.0)-p[0]*p[1]*xB[1];
  xBdot[1] = xB[0]+p[0]*xB[1];
}

