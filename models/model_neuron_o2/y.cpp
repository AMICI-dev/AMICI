
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void y_model_neuron_o2(double *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w) {
  y[0] = x[0];
  y[1] = x[2];
  y[2] = x[4];
  y[3] = x[6];
  y[4] = x[8];
}

} // namespace model_model_neuron_o2

} // namespace amici

