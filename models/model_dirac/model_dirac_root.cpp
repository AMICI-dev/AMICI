
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_dirac{

void root_model_dirac(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  root[0] = -t+p[1];
  root[1] = t-p[1];
}

} // namespace model_model_dirac

} // namespace amici

