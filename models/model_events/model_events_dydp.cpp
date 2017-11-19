
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dydp_model_events(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip) {
switch (ip) {
  case 3: {
  dydp[0] = x[0]+x[1]+x[2];

  } break;

}
}

