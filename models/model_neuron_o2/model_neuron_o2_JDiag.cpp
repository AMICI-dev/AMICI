
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void JDiag_model_neuron_o2(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JDiag[0+0*10] = x[0]*(2.0/2.5E1)+5.0;
  JDiag[1+0*10] = -p[0];
  JDiag[2+0*10] = w[1];
  JDiag[3+0*10] = -p[0];
  JDiag[4+0*10] = w[1];
  JDiag[5+0*10] = -p[0];
  JDiag[6+0*10] = w[1];
  JDiag[7+0*10] = -p[0];
  JDiag[8+0*10] = w[1];
  JDiag[9+0*10] = -p[0];
}

} // namespace model_model_neuron_o2

} // namespace amici

