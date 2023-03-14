
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void dJydsigma_model_neuron_o2(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  dJydsigma[0+0*5] = 1.0/(sigmay[0]*sigmay[0]*sigmay[0])*pow(my[0]-y[0],2.0)*-1.0+1.0/sigmay[0];
  dJydsigma[1+0*5] = y[1]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[2+0*5] = y[2]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[3+0*5] = y[3]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
  dJydsigma[4+0*5] = y[4]*1.0/(sigmay[0]*sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*1.0;
    break;
}
}

} // namespace model_model_neuron_o2

} // namespace amici

