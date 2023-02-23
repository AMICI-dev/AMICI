
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void root_model_calvetti(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx) {
  root[0] = -t+1.0E1;
  root[1] = -t+1.2E1;
  root[2] = t-1.0E1;
  root[3] = t-1.2E1;
}

} // namespace model_model_calvetti

} // namespace amici

