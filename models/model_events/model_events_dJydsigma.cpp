
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void dJydsigma_model_events(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  dJydsigma[0] = 1.0/(sigmay[0]*sigmay[0]*sigmay[0])*pow(my[0]-y[0]*1.0,2.0)*-1.0+1.0/sigmay[0];
    break;
}
}

} // namespace model_model_events

} // namespace amici

