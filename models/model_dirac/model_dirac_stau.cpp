
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void stau_model_dirac(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) {
switch (ip) {
  case 1: {
    switch(ie) { 
        case 0: {
  stau[0] = 1.0;

        } break;

        case 1: {
  stau[0] = 1.0;

        } break;

    } 

  } break;

}
}

