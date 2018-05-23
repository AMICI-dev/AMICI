
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void sxdot_model_jakstat_adjoint(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp) {
  sxdot[0] = dxdotdp[0]+J[0]*sx[0]+J[16]*sx[8];
  sxdot[1] = dxdotdp[1]+J[1]*sx[0]+J[2]*sx[1];
  sxdot[2] = dxdotdp[2]+J[3]*sx[1]+J[4]*sx[2];
  sxdot[3] = dxdotdp[3]+J[5]*sx[2]+J[6]*sx[3];
  sxdot[4] = dxdotdp[4]+J[7]*sx[3]+J[8]*sx[4];
  sxdot[5] = dxdotdp[5]+J[9]*sx[4]+J[10]*sx[5];
  sxdot[6] = dxdotdp[6]+J[11]*sx[5]+J[12]*sx[6];
  sxdot[7] = dxdotdp[7]+J[13]*sx[6]+J[14]*sx[7];
  sxdot[8] = dxdotdp[8]+J[15]*sx[7]+J[17]*sx[8];
}

