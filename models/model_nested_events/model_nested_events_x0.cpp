
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void x0_model_nested_events(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = p[0];
}

