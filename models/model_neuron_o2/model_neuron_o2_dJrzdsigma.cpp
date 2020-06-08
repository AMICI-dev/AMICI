
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void dJrzdsigma_model_neuron_o2(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) {
switch(iz){
    case 0:
  dJrzdsigma[0+0*5] = (rz[0]*rz[0])*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*-1.0+1.0/sigmaz[0];
  dJrzdsigma[1+0*5] = rz[0]*rz[1]*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*-2.0;
  dJrzdsigma[2+0*5] = rz[0]*rz[2]*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*-2.0;
  dJrzdsigma[3+0*5] = rz[0]*rz[3]*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*-2.0;
  dJrzdsigma[4+0*5] = rz[0]*rz[4]*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*-2.0;
    break;
}
}

} // namespace model_model_neuron_o2

} // namespace amici

