
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void drzdx_model_events(double *drzdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  drzdx[0+0*1] = 1.0;
  drzdx[0+2*1] = -1.0;
}

