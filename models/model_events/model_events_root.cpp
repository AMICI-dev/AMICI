
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void root_model_events(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  root[0] = x[1]-x[2];
  root[1] = x[0]-x[2];
  root[2] = t-4.0;
  root[3] = t-p[3];
}

