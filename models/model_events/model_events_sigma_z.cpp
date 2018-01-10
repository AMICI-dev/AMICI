
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void sigma_z_model_events(double *sigmaz, const realtype t, const realtype *p, const realtype *k) {
  sigmaz[0] = 1.0;
  sigmaz[1] = 1.0;
}

