
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void z_model_neuron(double *z, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
    switch(ie) { 
        case 0: {
  z[0] = t;

        } break;

    } 
}

