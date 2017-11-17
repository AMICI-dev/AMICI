
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dJzdsigma_model_neuron(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
switch(iz){
    case 0:
  dJzdsigma[0] = 1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*pow(mz[0]-z[0]*1.0,2.0)*-1.0+1.0/sigmaz[0];
    break;
}
}

