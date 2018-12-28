
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void dwdx_model_neuron_o2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl) {
  dwdx[0] = 2.0/2.5E1;
  dwdx[1] = dwdx[0];
}

