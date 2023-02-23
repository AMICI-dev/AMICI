
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_robertson{

void sigmay_model_robertson(double *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y) {
  sigmay[0] = 1.0;
  sigmay[1] = 1.0;
  sigmay[2] = 1.0;
}

} // namespace model_model_robertson

} // namespace amici

