
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void stau_model_events(double *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *sx, const int ip, const int ie) {
switch (ip) {
  case 0: {
    switch(ie) { 
        case 0: {
  tdata->stau[ip] = (sx_tmp[1]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  tdata->stau[ip] = (sx_tmp[0]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);

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
  tdata->stau[ip] = (sx_tmp[1]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  tdata->stau[ip] = (sx_tmp[0]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);

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
  tdata->stau[ip] = (sx_tmp[1]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  tdata->stau[ip] = (sx_tmp[0]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);

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
  tdata->stau[ip] = (sx_tmp[1]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->p[2]*x_tmp[1]-tdata->p[1]*x_tmp[0]*exp(t*(-1.0/1.0E1)));

        } break;

        case 1: {
  tdata->stau[ip] = (sx_tmp[0]-sx_tmp[2])/(tdata->h[2]-x_tmp[2]+tdata->h[3]*tdata->p[0]*x_tmp[0]);

        } break;

        case 2: {

        } break;

        case 3: {
  tdata->stau[ip] = 1.0;

        } break;

    } 

  } break;

}
}

