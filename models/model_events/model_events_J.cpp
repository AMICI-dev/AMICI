
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void J_model_events(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0+0*3] = -h[3]*p[0];
  J[1+0*3] = p[1]*exp(t*(-1.0/1.0E1));
  J[1+1*3] = -p[2];
  J[2+2*3] = -1.0;
}

