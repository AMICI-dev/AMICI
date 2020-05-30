
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void dJydy_model_calvetti(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
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
    case 3:
  dJydy[0+3*1] = 1.0/(sigmay[3]*sigmay[3])*(my[3]*2.0-y[3]*2.0)*-5.0E-1;
    break;
    case 4:
  dJydy[0+4*1] = 1.0/(sigmay[4]*sigmay[4])*(my[4]*2.0-y[4]*2.0)*-5.0E-1;
    break;
    case 5:
  dJydy[0+5*1] = 1.0/(sigmay[5]*sigmay[5])*(my[5]*2.0-y[5]*2.0)*-5.0E-1;
    break;
}
}

} // namespace model_model_calvetti

} // namespace amici

