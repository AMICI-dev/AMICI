
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void Jz_model_events(double *nllh, const int iz, const realtype *p, const realtype *k, const double *z, const double *sigmaz, const double *mz) {
switch(iz){
    case 0:
  nllh[0] = amici::log((sigmaz[0]*sigmaz[0])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmaz[0]*sigmaz[0])*pow(mz[0]-z[0],2.0)*5.0E-1;
    break;
    case 1:
  nllh[0] = amici::log((sigmaz[1]*sigmaz[1])*3.141592653589793*2.0)*5.0E-1+1.0/(sigmaz[1]*sigmaz[1])*pow(mz[1]-z[1],2.0)*5.0E-1;
    break;
}
}

} // namespace model_model_events

} // namespace amici

