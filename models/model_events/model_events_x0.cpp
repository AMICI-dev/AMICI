
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void x0_model_events(realtype *x0, const realtype t, const realtype *p, const realtype *k) {
  x0[0] = k[0];
  x0[1] = k[1];
  x0[2] = k[2];
}

