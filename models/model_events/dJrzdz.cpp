
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void dJrzdz_model_events(double *dJrzdz, const int iz, const realtype *p, const realtype *k, const double *rz, const double *sigmaz) {
switch(iz){
    case 0:
  dJrzdz[0+0*1] = rz[0]*1.0/(sigmaz[0]*sigmaz[0])*1.0;
    break;
    case 1:
  dJrzdz[0+1*1] = rz[1]*1.0/(sigmaz[1]*sigmaz[1])*1.0;
    break;
}
}

} // namespace model_model_events

} // namespace amici

