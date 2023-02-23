
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron{

void sx0_model_neuron(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) {
switch (ip) {
  case 1: {
  sx0[1] = k[0];

  } break;

}
}

} // namespace model_model_neuron

} // namespace amici

