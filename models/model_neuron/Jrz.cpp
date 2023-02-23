
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron{

void Jrz_model_neuron(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) {
switch(iz){
    case 0:
  nllh[0] = amici::log((sigmaz[0]*sigmaz[0])*3.141592653589793*2.0)*5.0E-1+(rz[0]*rz[0])*1.0/(sigmaz[0]*sigmaz[0])*5.0E-1;
    break;
}
}

} // namespace model_model_neuron

} // namespace amici

