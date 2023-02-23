
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void rz_model_events(double *rz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
    switch(ie) { 
        case 0: {
  rz[0] = x[1]-x[2];

        } break;

        case 1: {
  rz[1] = x[0]-x[2];

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

