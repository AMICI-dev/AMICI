
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron{

void JB_model_neuron(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0+0*2] = x[0]*(-2.0/2.5E1)-5.0;
  JB[0+1*2] = -p[0]*p[1];
  JB[1+0*2] = 1.0;
  JB[1+1*2] = p[0];
}

} // namespace model_model_neuron

} // namespace amici

