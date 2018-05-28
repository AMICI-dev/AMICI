
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void JvB_model_neuron(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx) {
  JvB[0] = -vB[0]*(x[0]*(2.0/2.5E1)+5.0)-p[0]*p[1]*vB[1];
  JvB[1] = vB[0]+p[0]*vB[1];
}

