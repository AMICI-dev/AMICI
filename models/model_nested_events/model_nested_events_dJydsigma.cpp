
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void dJydsigma_model_nested_events(double *dJydsigma, const int iy, const realtype *p, const realtype *k, const double *y, const double *sigmay, const double *my) {
switch(iy){
    case 0:
  dJydsigma[0] = 1.0/(sigmay[0]*sigmay[0]*sigmay[0])*pow(my[0]-y[0]*1.0,2.0)*-1.0+1.0/sigmay[0];
    break;
}
}

