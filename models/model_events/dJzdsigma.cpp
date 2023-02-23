
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void dJzdsigma_model_events(double *dJzdsigma, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
switch(iz){
    case 0:
  dJzdsigma[0+0*1] = 1.0/(sigmaz[0]*sigmaz[0]*sigmaz[0])*pow(mz[0]-z[0],2.0)*-1.0+1.0/sigmaz[0];
    break;
    case 1:
  dJzdsigma[0+1*1] = 1.0/(sigmaz[1]*sigmaz[1]*sigmaz[1])*pow(mz[1]-z[1],2.0)*-1.0+1.0/sigmaz[1];
    break;
}
}

} // namespace model_model_events

} // namespace amici

