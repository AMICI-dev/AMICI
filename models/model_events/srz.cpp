
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_events{

void srz_model_events(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[1]-sx[2];

        } break;

        case 1: {
  srz[1] = sx[0]-sx[2];

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

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[1]-sx[2];

        } break;

        case 1: {
  srz[1] = sx[0]-sx[2];

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

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[1]-sx[2];

        } break;

        case 1: {
  srz[1] = sx[0]-sx[2];

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

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[1]-sx[2];

        } break;

        case 1: {
  srz[1] = sx[0]-sx[2];

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

  } break;

}
}

} // namespace model_model_events

} // namespace amici

