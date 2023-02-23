
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void dJzdsigma_model_neuron_o2(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
switch(iz){
    case 0:
  dJzdsigma[0+0*5] = 1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*pow(mz[0]-z[0],2.0)*-1.0+1.0/sigmaz[0];
  dJzdsigma[1+0*5] = z[1]*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*(mz[0]*2.0-z[0]*2.0)*1.0;
  dJzdsigma[2+0*5] = z[2]*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*(mz[0]*2.0-z[0]*2.0)*1.0;
  dJzdsigma[3+0*5] = z[3]*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*(mz[0]*2.0-z[0]*2.0)*1.0;
  dJzdsigma[4+0*5] = z[4]*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*(mz[0]*2.0-z[0]*2.0)*1.0;
    break;
}
}

} // namespace model_model_neuron_o2

} // namespace amici

