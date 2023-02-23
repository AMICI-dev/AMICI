
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_robertson{

void x0_model_robertson(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = k[0];
}

} // namespace model_model_robertson

} // namespace amici

