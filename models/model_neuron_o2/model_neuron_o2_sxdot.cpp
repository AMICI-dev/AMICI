
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void sxdot_model_neuron_o2(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp) {
  sxdot[0] = dxdotdp[0]+J[0]*sx[0]+J[8]*sx[1];
  sxdot[1] = dxdotdp[1]+J[1]*sx[0]+J[9]*sx[1];
  sxdot[2] = dxdotdp[2]+J[2]*sx[0]+J[11]*sx[2]+J[13]*sx[3];
  sxdot[3] = dxdotdp[3]+J[3]*sx[0]+J[10]*sx[1]+J[12]*sx[2]+J[14]*sx[3];
  sxdot[4] = dxdotdp[4]+J[4]*sx[0]+J[15]*sx[4]+J[17]*sx[5];
  sxdot[5] = dxdotdp[5]+J[5]*sx[0]+J[16]*sx[4]+J[18]*sx[5];
  sxdot[6] = dxdotdp[6]+J[6]*sx[0]+J[19]*sx[6]+J[21]*sx[7];
  sxdot[7] = dxdotdp[7]+J[20]*sx[6]+J[22]*sx[7];
  sxdot[8] = dxdotdp[8]+J[7]*sx[0]+J[23]*sx[8]+J[25]*sx[9];
  sxdot[9] = dxdotdp[9]+J[24]*sx[8]+J[26]*sx[9];
}

