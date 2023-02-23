
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void z_model_neuron_o2(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
    switch(ie) { 
        case 0: {
  z[0] = t;
  z[1] = -x[2]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  z[2] = -x[4]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  z[3] = -x[6]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  z[4] = -x[8]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 
}

} // namespace model_model_neuron_o2

} // namespace amici

