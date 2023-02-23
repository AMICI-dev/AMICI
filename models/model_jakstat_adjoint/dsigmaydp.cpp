
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void dsigmaydp_model_jakstat_adjoint(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip) {
switch (ip) {
  case 14: {
  dsigmaydp[0] = 1.0;

  } break;

  case 15: {
  dsigmaydp[1] = 1.0;

  } break;

  case 16: {
  dsigmaydp[2] = 1.0;

  } break;

}
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

