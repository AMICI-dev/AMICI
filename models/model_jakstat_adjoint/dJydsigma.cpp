
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void dJydsigma_model_jakstat_adjoint(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
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
}
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

