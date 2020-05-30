
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_robertson{

void Jy_model_robertson(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
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
}
}

} // namespace model_model_robertson

} // namespace amici

