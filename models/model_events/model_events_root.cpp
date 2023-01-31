
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void root_model_events(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) {
  root[0] = x[1]-x[2];
  root[1] = x[0]-x[2];
  root[2] = -t+4.0;
  root[3] = -t+p[3];
  root[4] = t-4.0;
  root[5] = t-p[3];
}

} // namespace model_model_events

} // namespace amici

