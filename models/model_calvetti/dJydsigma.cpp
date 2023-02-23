
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void dJydsigma_model_calvetti(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  dJydsigma[0+0*1] = 1.0/(sigmay[0]*sigmay[0]*sigmay[0])*pow(my[0]-y[0],2.0)*-1.0+1.0/sigmay[0];
    break;
    case 1:
  dJydsigma[0+1*1] = 1.0/(sigmay[1]*sigmay[1]*sigmay[1])*pow(my[1]-y[1],2.0)*-1.0+1.0/sigmay[1];
    break;
    case 2:
  dJydsigma[0+2*1] = 1.0/(sigmay[2]*sigmay[2]*sigmay[2])*pow(my[2]-y[2],2.0)*-1.0+1.0/sigmay[2];
    break;
    case 3:
  dJydsigma[0+3*1] = 1.0/(sigmay[3]*sigmay[3]*sigmay[3])*pow(my[3]-y[3],2.0)*-1.0+1.0/sigmay[3];
    break;
    case 4:
  dJydsigma[0+4*1] = 1.0/(sigmay[4]*sigmay[4]*sigmay[4])*pow(my[4]-y[4],2.0)*-1.0+1.0/sigmay[4];
    break;
    case 5:
  dJydsigma[0+5*1] = 1.0/(sigmay[5]*sigmay[5]*sigmay[5])*pow(my[5]-y[5],2.0)*-1.0+1.0/sigmay[5];
    break;
}
}

} // namespace model_model_calvetti

} // namespace amici

