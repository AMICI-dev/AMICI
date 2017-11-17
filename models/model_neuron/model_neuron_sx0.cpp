
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void sx0_model_neuron(realtype *sx0, const realtype t,const realtype *x0, const realtype *p, const realtype *k, const int ip) {
switch (ip) {
  case 1: {
  sx0[1] = k[0];

  } break;

}
}

