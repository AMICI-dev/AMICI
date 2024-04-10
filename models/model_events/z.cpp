
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void z_model_events(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
    switch(ie) { 
        case 0: {
  z[0] = t;

        } break;

        case 1: {
  z[1] = t;

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

        case 4: {

        } break;

        case 5: {

        } break;

    } 
}

} // namespace model_model_events

} // namespace amici

