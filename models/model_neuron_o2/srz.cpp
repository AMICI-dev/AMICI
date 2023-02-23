
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void srz_model_neuron_o2(double *srz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[0];
  srz[1] = sx[2];
  srz[2] = sx[4];
  srz[3] = sx[6];
  srz[4] = sx[8];

        } break;

    } 

  } break;

  case 1: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[0];
  srz[1] = sx[2];
  srz[2] = sx[4];
  srz[3] = sx[6];
  srz[4] = sx[8];

        } break;

    } 

  } break;

  case 2: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[0];
  srz[1] = sx[2];
  srz[2] = sx[4];
  srz[3] = sx[6];
  srz[4] = sx[8];

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  srz[0] = sx[0];
  srz[1] = sx[2];
  srz[2] = sx[4];
  srz[3] = sx[6];
  srz[4] = sx[8];

        } break;

    } 

  } break;

}
}

} // namespace model_model_neuron_o2

} // namespace amici

