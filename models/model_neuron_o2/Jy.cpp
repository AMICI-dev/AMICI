
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void Jy_model_neuron_o2(double *nllh, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  nllh[0] = amici::log((sigmay[0]*sigmay[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmay[0]*sigmay[0])*pow(my[0]-y[0],2.0)*5.0E-1;
  nllh[1] = y[1]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[2] = y[2]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[3] = y[3]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
  nllh[4] = y[4]*1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
    break;
}
}

} // namespace model_model_neuron_o2

} // namespace amici

