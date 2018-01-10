
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void drzdx_model_events(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  drzdx[0+0*1] = 1.0;
  drzdx[0+2*1] = -1.0;
}

