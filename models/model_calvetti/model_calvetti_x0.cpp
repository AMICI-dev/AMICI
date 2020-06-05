
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void x0_model_calvetti(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = k[0];
  x0[1] = k[2];
  x0[2] = k[4];
  x0[3] = 1.0;
  x0[4] = 1.0;
  x0[5] = 1.0;
}

} // namespace model_model_calvetti

} // namespace amici

