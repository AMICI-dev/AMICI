
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void drzdx_model_neuron_o2(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  drzdx[0+0*5] = 1.0;
  drzdx[1+2*5] = 1.0;
  drzdx[2+4*5] = 1.0;
  drzdx[3+6*5] = 1.0;
  drzdx[4+8*5] = 1.0;
}

} // namespace model_model_neuron_o2

} // namespace amici

