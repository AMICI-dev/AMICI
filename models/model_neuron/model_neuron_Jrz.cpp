
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void Jrz_model_neuron(double *nllh, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) {
switch(iz){
    case 0:
  nllh[0] = amici::amilog((sigmaz[0]*sigmaz[0])*3.141592653589793*2.0)*5.0E-1+(rz[0]*rz[0])*1.0/(sigmaz[0]*sigmaz[0])*5.0E-1;
    break;
}
}

