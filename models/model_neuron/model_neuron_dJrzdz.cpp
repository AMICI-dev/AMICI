
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dJrzdz_model_neuron(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) {
switch(iz){
    case 0:
  dJrzdz[0] = rz[0]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
    break;
}
}

