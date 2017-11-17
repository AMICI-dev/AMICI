
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void Jz_model_neuron(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
switch(iz){
    case 0:
  nllh[0] = amici::amilog((sigmaz[0]*sigmaz[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmaz[0]*sigmaz[0])*pow(mz[0]-z[0],2.0)*5.0E-1;
    break;
}
}

