
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void dzdx_model_neuron_o2(double *dzdx, const int ie, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  dzdx[0+0*5] = -1.0/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  dzdx[1+0*5] = x[2]*(x[0]*(2.0/2.5E1)+5.0)*1.0/pow(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2,2.0);
  dzdx[1+1*5] = -x[2]*1.0/pow(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2,2.0);
  dzdx[1+2*5] = -1.0/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  dzdx[2+0*5] = x[4]*(x[0]*(2.0/2.5E1)+5.0)*1.0/pow(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2,2.0);
  dzdx[2+1*5] = -x[4]*1.0/pow(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2,2.0);
  dzdx[2+4*5] = -1.0/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  dzdx[3+0*5] = x[6]*(x[0]*(2.0/2.5E1)+5.0)*1.0/pow(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2,2.0);
  dzdx[3+1*5] = -x[6]*1.0/pow(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2,2.0);
  dzdx[3+6*5] = -1.0/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  dzdx[4+0*5] = x[8]*(x[0]*(2.0/2.5E1)+5.0)*1.0/pow(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2,2.0);
  dzdx[4+1*5] = -x[8]*1.0/pow(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2,2.0);
  dzdx[4+8*5] = -1.0/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
}

} // namespace model_model_neuron_o2

} // namespace amici

