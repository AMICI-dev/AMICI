
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_nested_events{

void deltasx_model_nested_events(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau) {
switch (ip) {
  case 0: {
              switch(ie) { 
              case 0: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 1: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 2: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 3: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              } 

  } break;

  case 1: {
              switch(ie) { 
              case 0: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 1: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 2: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0])+1.0;

              } break;

              case 3: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0])+1.0;

              } break;

              } 

  } break;

  case 2: {
              switch(ie) { 
              case 0: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 1: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 2: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 3: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              } 

  } break;

  case 3: {
              switch(ie) { 
              case 0: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 1: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 2: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 3: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              } 

  } break;

  case 4: {
              switch(ie) { 
              case 0: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 1: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 2: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              case 3: {
  deltasx[0] = -stau[0]*(xdot[0]-xdot_old[0]);

              } break;

              } 

  } break;

}
}

} // namespace model_model_nested_events

} // namespace amici

