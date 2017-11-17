
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dsigma_ydp_model_jakstat_adjoint_o2(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) {
switch (ip) {
  case 14: {
  dsigmaydp[ip*54 + 0] = 1.0;

  } break;

  case 15: {
  dsigmaydp[ip*54 + 1] = 1.0;

  } break;

  case 16: {
  dsigmaydp[ip*54 + 2] = 1.0;

  } break;

}
}

