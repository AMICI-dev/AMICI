
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void sxdot_model_neuron_o2(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp) {
  sxdot[0] = dxdotdp[0]+J[0]*sx[0]+J[10]*sx[1];
  sxdot[1] = dxdotdp[1]+J[1]*sx[0]+J[11]*sx[1];
  sxdot[2] = dxdotdp[2]+J[2]*sx[0]+J[22]*sx[2]+J[32]*sx[3];
  sxdot[3] = dxdotdp[3]+J[3]*sx[0]+J[13]*sx[1]+J[23]*sx[2]+J[33]*sx[3];
  sxdot[4] = dxdotdp[4]+J[4]*sx[0]+J[44]*sx[4]+J[54]*sx[5];
  sxdot[5] = dxdotdp[5]+J[5]*sx[0]+J[45]*sx[4]+J[55]*sx[5];
  sxdot[6] = dxdotdp[6]+J[6]*sx[0]+J[66]*sx[6]+J[76]*sx[7];
  sxdot[7] = dxdotdp[7]+J[67]*sx[6]+J[77]*sx[7];
  sxdot[8] = dxdotdp[8]+J[8]*sx[0]+J[88]*sx[8]+J[98]*sx[9];
  sxdot[9] = dxdotdp[9]+J[89]*sx[8]+J[99]*sx[9];
}

