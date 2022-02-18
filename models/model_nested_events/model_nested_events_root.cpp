
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_nested_events{

void root_model_nested_events(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl) {
  root[0] = -x[0]+1.0;
  root[1] = x[0]-1.0;
  root[2] = t-p[2];
  root[3] = -t+p[2];
}

} // namespace model_model_nested_events

} // namespace amici

