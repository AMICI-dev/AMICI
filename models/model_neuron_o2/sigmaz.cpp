
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void sigmaz_model_neuron_o2(double *sigmaz, const realtype t, const realtype *p, const realtype *k) {
  sigmaz[0] = 1.0;
}

} // namespace model_model_neuron_o2

} // namespace amici

