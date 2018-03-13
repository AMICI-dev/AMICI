
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dJrzdsigma_model_neuron(double *dJrzdsigma, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) {
switch(iz){
    case 0:
  dJrzdsigma[0] = (rz[0]*rz[0])*1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*-1.0+1.0/sigmaz[0];
    break;
}
}

