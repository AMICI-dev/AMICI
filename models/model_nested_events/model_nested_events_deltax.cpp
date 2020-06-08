
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_nested_events{

void deltax_model_nested_events(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) {
              switch(ie) { 
              case 2: {
  deltax[0] = p[1];

              } break;

              case 3: {
  deltax[0] = p[1];

              } break;

              } 
}

} // namespace model_model_nested_events

} // namespace amici

