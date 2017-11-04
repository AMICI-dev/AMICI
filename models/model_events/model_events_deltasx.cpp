
#include <include/symbolic_functions.h>
#include "model_events_w.h"

using namespace model_events;

void deltasx_model_events(double *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau) {
switch (udata->plist[ip]) {
  case 0: {
              switch(ie) { 
              case 2: {
  tdata->deltasx[ip*model->nx + 2] = -tdata->stau[ip]*(xdot_tmp[2]-xdot_old_tmp[2]);

              } break;

              case 3: {
  tdata->deltasx[ip*model->nx + 0] = -tdata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);

              } break;

              } 

  } break;

  case 1: {
              switch(ie) { 
              case 2: {
  tdata->deltasx[ip*model->nx + 2] = -tdata->stau[ip]*(xdot_tmp[2]-xdot_old_tmp[2]);

              } break;

              case 3: {
  tdata->deltasx[ip*model->nx + 0] = -tdata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);

              } break;

              } 

  } break;

  case 2: {
              switch(ie) { 
              case 2: {
  tdata->deltasx[ip*model->nx + 2] = -tdata->stau[ip]*(xdot_tmp[2]-xdot_old_tmp[2]);

              } break;

              case 3: {
  tdata->deltasx[ip*model->nx + 0] = -tdata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);

              } break;

              } 

  } break;

  case 3: {
              switch(ie) { 
              case 2: {
  tdata->deltasx[ip*model->nx + 2] = -tdata->stau[ip]*(xdot_tmp[2]-xdot_old_tmp[2]);

              } break;

              case 3: {
  tdata->deltasx[ip*model->nx + 0] = -tdata->stau[ip]*(xdot_tmp[0]-xdot_old_tmp[0]);

              } break;

              } 

  } break;

}
}

