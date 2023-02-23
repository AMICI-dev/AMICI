
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron{

void dJzdsigma_model_neuron(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
switch(iz){
    case 0:
  dJzdsigma[0] = 1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*pow(mz[0]-z[0]*1.0,2.0)*-1.0+1.0/sigmaz[0];
    break;
}
}

} // namespace model_model_neuron

} // namespace amici

