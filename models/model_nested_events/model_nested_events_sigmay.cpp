
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void sigmay_model_nested_events(double *sigmay, const realtype t, const realtype *p, const realtype *k) {
  sigmay[0] = 1.0;
}

