
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dsigma_ydp_model_jakstat_adjoint(double *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const int ip) {
switch (ip) {
  case 14: {
  dsigmaydp[ip*3 + 0] = 1.0;

  } break;

  case 15: {
  dsigmaydp[ip*3 + 1] = 1.0;

  } break;

  case 16: {
  dsigmaydp[ip*3 + 2] = 1.0;

  } break;

}
}

