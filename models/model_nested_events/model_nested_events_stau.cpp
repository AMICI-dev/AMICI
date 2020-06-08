
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_nested_events{

void stau_model_nested_events(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 1: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 1: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 1: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 2: {
  stau[0] = 1.0;

        } break;

        case 3: {
  stau[0] = 1.0;

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 1: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 4: {
    switch(ie) { 
        case 0: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 1: {
  stau[0] = sx[0]/(p[4]*x[0]+p[3]*x[0]*(h[0]-1.0));

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

}
}

} // namespace model_model_nested_events

} // namespace amici

