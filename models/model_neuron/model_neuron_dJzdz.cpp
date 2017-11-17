
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dJzdz_model_neuron(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
switch(iz){
    case 0:
  dJzdz[0] = 1.0/(sigmaz[0]*sigmaz[0])*(mz[0]*2.0-z[0]*2.0)*-5.0E-1;
    break;
}
}

