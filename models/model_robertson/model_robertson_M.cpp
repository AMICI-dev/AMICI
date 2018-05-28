
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void M_model_robertson(realtype *M, const realtype t, const realtype *x, const realtype *p, const realtype *k) {
  M[0+0*3] = 1.0;
  M[1+1*3] = 1.0;
}

