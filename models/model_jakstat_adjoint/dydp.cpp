
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_jakstat_adjoint{

void dydp_model_jakstat_adjoint(double *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 4: {
  dydp[0] = -1.0/(p[4]*p[4])*p[13]*(x[1]+x[2]*2.0);
  dydp[1] = -1.0/(p[4]*p[4])*p[12]*(x[0]+x[1]+x[2]*2.0);

  } break;

  case 5: {
  dydp[2] = amici::Dspline_pos(4,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);

  } break;

  case 6: {
  dydp[2] = amici::Dspline_pos(6,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);

  } break;

  case 7: {
  dydp[2] = amici::Dspline_pos(8,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);

  } break;

  case 8: {
  dydp[2] = amici::Dspline_pos(10,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);

  } break;

  case 9: {
  dydp[2] = amici::Dspline_pos(12,t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);

  } break;

  case 10: {
  dydp[1] = 1.0;

  } break;

  case 11: {
  dydp[0] = 1.0;

  } break;

  case 12: {
  dydp[1] = (x[0]+x[1]+x[2]*2.0)/p[4];

  } break;

  case 13: {
  dydp[0] = (x[1]+x[2]*2.0)/p[4];

  } break;

}
}

} // namespace model_model_jakstat_adjoint

} // namespace amici

