
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void stau_model_neuron(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  stau[ip] = -sx[0]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  stau[ip] = -sx[0]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  stau[ip] = -sx[0]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  stau[ip] = -sx[0]/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);

        } break;

    } 

  } break;

}
}

