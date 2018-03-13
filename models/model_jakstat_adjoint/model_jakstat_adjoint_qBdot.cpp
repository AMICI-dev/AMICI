
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void qBdot_model_jakstat_adjoint(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  qBdot[0] = w[0]*x[0]*xB[0]-w[0]*x[0]*xB[1];

  } break;

  case 1: {
  qBdot[0] = w[1]*xB[1]*2.0-w[1]*xB[2];

  } break;

  case 2: {
  qBdot[0] = x[2]*xB[2]-(k[0]*x[2]*xB[3])/k[1];

  } break;

  case 3: {
  qBdot[0] = x[3]*xB[3]-xB[5]*(x[4]-x[5])-xB[6]*(x[5]-x[6])-xB[7]*(x[6]-x[7])-xB[8]*(x[7]-x[8])-xB[4]*(x[3]*2.0-x[4])-(k[1]*x[8]*xB[0])/k[0];

  } break;

  case 5: {
  qBdot[0] = p[0]*x[0]*xB[0]*dwdp[0]-p[0]*x[0]*xB[1]*dwdp[0];

  } break;

  case 6: {
  qBdot[0] = p[0]*x[0]*xB[0]*dwdp[1]-p[0]*x[0]*xB[1]*dwdp[1];

  } break;

  case 7: {
  qBdot[0] = p[0]*x[0]*xB[0]*dwdp[2]-p[0]*x[0]*xB[1]*dwdp[2];

  } break;

  case 8: {
  qBdot[0] = p[0]*x[0]*xB[0]*dwdp[3]-p[0]*x[0]*xB[1]*dwdp[3];

  } break;

  case 9: {
  qBdot[0] = p[0]*x[0]*xB[0]*dwdp[4]-p[0]*x[0]*xB[1]*dwdp[4];

  } break;

}
}

