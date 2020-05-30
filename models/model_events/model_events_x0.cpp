
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void x0_model_events(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = k[0];
  x0[1] = k[1];
  x0[2] = k[2];
}

} // namespace model_model_events

} // namespace amici

