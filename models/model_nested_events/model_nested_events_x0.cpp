
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_nested_events{

void x0_model_nested_events(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = p[0];
}

} // namespace model_model_nested_events

} // namespace amici

