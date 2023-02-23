
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void dJrzdz_model_neuron_o2(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) {
switch(iz){
    case 0:
  dJrzdz[0+0*5] = rz[0]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
  dJrzdz[1+0*5] = rz[1]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
  dJrzdz[1+1*5] = rz[0]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
  dJrzdz[2+0*5] = rz[2]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
  dJrzdz[2+2*5] = rz[0]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
  dJrzdz[3+0*5] = rz[3]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
  dJrzdz[3+3*5] = rz[0]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
  dJrzdz[4+0*5] = rz[4]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
  dJrzdz[4+4*5] = rz[0]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
    break;
}
}

} // namespace model_model_neuron_o2

} // namespace amici

