
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void dJzdz_model_events(double *dJzdz, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
switch(iz){
    case 0:
  dJzdz[0+0*1] = 1.0/(sigmaz[0]*sigmaz[0])*(mz[0]*2.0-z[0]*2.0)*-5.0E-1;
    break;
    case 1:
  dJzdz[0+1*1] = 1.0/(sigmaz[1]*sigmaz[1])*(mz[1]*2.0-z[1]*2.0)*-5.0E-1;
    break;
}
}

} // namespace model_model_events

} // namespace amici

