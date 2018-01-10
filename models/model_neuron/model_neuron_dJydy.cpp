
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dJydy_model_neuron(double *dJydy, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  dJydy[0] = 1.0/(sigmay[0]*sigmay[0])*(my[0]*2.0-y[0]*2.0)*-5.0E-1;
    break;
}
}

