
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void sz_model_events(double *sz, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  sz[0] = (sx[1]-sx[2])/(h[2]-x[2]+p[2]*x[1]-p[1]*x[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  sz[1] = (sx[0]-sx[2])/(h[2]-x[2]+h[3]*p[0]*x[0]);

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
  sz[0] = (sx[1]-sx[2])/(h[2]-x[2]+p[2]*x[1]-p[1]*x[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  sz[1] = (sx[0]-sx[2])/(h[2]-x[2]+h[3]*p[0]*x[0]);

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
  sz[0] = (sx[1]-sx[2])/(h[2]-x[2]+p[2]*x[1]-p[1]*x[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  sz[1] = (sx[0]-sx[2])/(h[2]-x[2]+h[3]*p[0]*x[0]);

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

  case 3: {
    switch(ie) { 
        case 0: {
  sz[0] = (sx[1]-sx[2])/(h[2]-x[2]+p[2]*x[1]-p[1]*x[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  sz[1] = (sx[0]-sx[2])/(h[2]-x[2]+h[3]*p[0]*x[0]);

        } break;

        case 2: {

        } break;

        case 3: {

        } break;

    } 

  } break;

}
}

