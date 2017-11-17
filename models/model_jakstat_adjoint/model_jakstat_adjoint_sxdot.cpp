
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void sxdot_model_jakstat_adjoint(realtype *sxdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *sx, const realtype *w, const realtype *dwdx, const realtype *J, const realtype *dxdotdp) {
  sxdot[0] = dxdotdp[0]+J[0]*sx[0]+J[72]*sx[8];
  sxdot[1] = dxdotdp[1]+J[1]*sx[0]+J[10]*sx[1];
  sxdot[2] = dxdotdp[2]+J[11]*sx[1]+J[20]*sx[2];
  sxdot[3] = dxdotdp[3]+J[21]*sx[2]+J[30]*sx[3];
  sxdot[4] = dxdotdp[4]+J[31]*sx[3]+J[40]*sx[4];
  sxdot[5] = dxdotdp[5]+J[41]*sx[4]+J[50]*sx[5];
  sxdot[6] = dxdotdp[6]+J[51]*sx[5]+J[60]*sx[6];
  sxdot[7] = dxdotdp[7]+J[61]*sx[6]+J[70]*sx[7];
  sxdot[8] = dxdotdp[8]+J[71]*sx[7]+J[80]*sx[8];
}

