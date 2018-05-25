
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void w_model_robertson(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  w[0] = p[1]*x[1]*x[2];
}

