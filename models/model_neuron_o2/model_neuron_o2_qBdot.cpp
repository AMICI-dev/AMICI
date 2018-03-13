
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void qBdot_model_neuron_o2(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  qBdot[0] = xB[1]*(x[1]-p[1]*x[0]);
  qBdot[1] = xB[3]*(x[1]-p[1]*x[0])+xB[1]*(x[3]-p[1]*x[2]);
  qBdot[2] = -xB[1]*(x[0]-x[5]+p[1]*x[4])+xB[5]*(x[1]-p[1]*x[0]);
  qBdot[3] = xB[7]*(x[1]-p[1]*x[0])+xB[1]*(x[7]-p[1]*x[6]);
  qBdot[4] = xB[9]*(x[1]-p[1]*x[0])+xB[1]*(x[9]-p[1]*x[8]);

  } break;

  case 1: {
  qBdot[0] = -p[0]*x[0]*xB[1];
  qBdot[1] = -xB[1]*(x[0]+p[0]*x[2])-p[0]*x[0]*xB[3];
  qBdot[2] = -p[0]*x[0]*xB[5]-p[0]*x[4]*xB[1];
  qBdot[3] = -p[0]*x[0]*xB[7]-p[0]*x[6]*xB[1];
  qBdot[4] = -p[0]*x[0]*xB[9]-p[0]*x[8]*xB[1];

  } break;

}
}

