
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void deltaqB_model_neuron(double *deltaqB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *xB) {
switch (ip) {
  case 2: {
              switch(ie) { 
              case 0: {
  deltaqB[ip+0] = xB[0];

              } break;

              } 

  } break;

  case 3: {
              switch(ie) { 
              case 0: {
  deltaqB[ip+0] = -xB[1];

              } break;

              } 

  } break;

}
}

