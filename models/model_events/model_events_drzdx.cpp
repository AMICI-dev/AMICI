
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void drzdx_model_events(double *drzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  drzdx[0+1*2] = 1.0;
  drzdx[0+2*2] = -1.0;
  drzdx[1+0*2] = 1.0;
  drzdx[1+2*2] = -1.0;
}

} // namespace model_model_events

} // namespace amici

