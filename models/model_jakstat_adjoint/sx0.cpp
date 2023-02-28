
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void sx0_model_jakstat_adjoint(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) {
switch (ip) {
  case 4: {
  sx0[0] = 1.0;

  } break;

}
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

