
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void sz_model_neuron(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  sz[0] = -sx[0]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  sz[0] = -sx[0]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  sz[0] = -sx[0]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  sz[0] = -sx[0]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

}
}

