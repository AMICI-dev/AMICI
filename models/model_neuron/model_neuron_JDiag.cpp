
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron{

void JDiag_model_neuron(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JDiag[0+0*2] = x[0]*(2.0/2.5E1)+5.0;
  JDiag[1+0*2] = -p[0];
}

} // namespace model_model_neuron

} // namespace amici

