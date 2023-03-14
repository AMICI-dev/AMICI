
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void x0_model_jakstat_adjoint(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = p[4];
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

