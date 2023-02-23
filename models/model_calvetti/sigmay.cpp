
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void sigmay_model_calvetti(double *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) {
  sigmay[0] = 1.0;
  sigmay[1] = 1.0;
  sigmay[2] = 1.0;
  sigmay[3] = 1.0;
  sigmay[4] = 1.0;
  sigmay[5] = 1.0;
}

} // namespace model_model_calvetti

} // namespace amici

