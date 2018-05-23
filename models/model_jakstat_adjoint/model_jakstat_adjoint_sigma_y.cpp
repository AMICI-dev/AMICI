
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void sigma_y_model_jakstat_adjoint(double *sigmay, const realtype t, const realtype *p, const realtype *k) {
  sigmay[0] = p[14];
  sigmay[1] = p[15];
  sigmay[2] = p[16];
}

