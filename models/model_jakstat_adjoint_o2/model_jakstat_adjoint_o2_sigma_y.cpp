
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void sigma_y_model_jakstat_adjoint_o2(double *sigmay, const realtype t, const realtype *p, const realtype *k) {
  sigmay[0] = p[14];
  sigmay[1] = p[15];
  sigmay[2] = p[16];
  sigmay[17] = 1.0;
  sigmay[35] = 1.0;
  sigmay[53] = 1.0;
}

