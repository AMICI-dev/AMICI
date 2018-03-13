
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dJydy_model_robertson(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  dJydy[0+0*1] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
    break;
    case 1:
  dJydy[0+1*1] = 1.0/(sigmay[1]*sigmay[1])*(my[1]*2.0-y[1]*2.0)*-5.0E-1;
    break;
    case 2:
  dJydy[0+2*1] = 1.0/(sigmay[2]*sigmay[2])*(my[2]*2.0-y[2]*2.0)*-5.0E-1;
    break;
}
}

