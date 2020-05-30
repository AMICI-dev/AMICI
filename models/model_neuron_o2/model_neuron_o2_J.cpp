
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void J_model_neuron_o2(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  J[0+0*10] = x[0]*(2.0/2.5E1)+5.0;
  J[0+1*10] = -1.0;
  J[1+0*10] = p[0]*p[1];
  J[1+1*10] = -p[0];
  J[2+0*10] = x[2]*dwdx[1];
  J[2+2*10] = w[1];
  J[2+3*10] = -1.0;
  J[3+0*10] = p[1];
  J[3+1*10] = -1.0;
  J[3+2*10] = p[0]*p[1];
  J[3+3*10] = -p[0];
  J[4+0*10] = x[4]*dwdx[1];
  J[4+4*10] = w[1];
  J[4+5*10] = -1.0;
  J[5+0*10] = p[0];
  J[5+4*10] = p[0]*p[1];
  J[5+5*10] = -p[0];
  J[6+0*10] = x[6]*dwdx[1];
  J[6+6*10] = w[1];
  J[6+7*10] = -1.0;
  J[7+6*10] = p[0]*p[1];
  J[7+7*10] = -p[0];
  J[8+0*10] = x[8]*dwdx[1];
  J[8+8*10] = w[1];
  J[8+9*10] = -1.0;
  J[9+8*10] = p[0]*p[1];
  J[9+9*10] = -p[0];
}

} // namespace model_model_neuron_o2

} // namespace amici

