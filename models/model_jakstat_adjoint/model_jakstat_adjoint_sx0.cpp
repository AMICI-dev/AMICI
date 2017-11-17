
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void sx0_model_jakstat_adjoint(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) {
switch (ip) {
  case 4: {
  sx0[0] = 1.0;

  } break;

}
}

