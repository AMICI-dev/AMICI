
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void Jy_model_calvetti(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  nllh[0] = amici::log((sigmay[0]*sigmay[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[0]*sigmay[0])*pow(my[0]-y[0],2.0)*5.0E-1;
    break;
    case 1:
  nllh[0] = amici::log((sigmay[1]*sigmay[1])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[1]*sigmay[1])*pow(my[1]-y[1],2.0)*5.0E-1;
    break;
    case 2:
  nllh[0] = amici::log((sigmay[2]*sigmay[2])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[2]*sigmay[2])*pow(my[2]-y[2],2.0)*5.0E-1;
    break;
    case 3:
  nllh[0] = amici::log((sigmay[3]*sigmay[3])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[3]*sigmay[3])*pow(my[3]-y[3],2.0)*5.0E-1;
    break;
    case 4:
  nllh[0] = amici::log((sigmay[4]*sigmay[4])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[4]*sigmay[4])*pow(my[4]-y[4],2.0)*5.0E-1;
    break;
    case 5:
  nllh[0] = amici::log((sigmay[5]*sigmay[5])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[5]*sigmay[5])*pow(my[5]-y[5],2.0)*5.0E-1;
    break;
}
}

} // namespace model_model_calvetti

} // namespace amici

