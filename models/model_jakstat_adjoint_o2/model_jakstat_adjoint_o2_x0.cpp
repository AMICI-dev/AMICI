
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint_o2{

void x0_model_jakstat_adjoint_o2(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = p[4];
  x0[45] = 1.0;
}

} // namespace model_model_jakstat_adjoint_o2

} // namespace amici

