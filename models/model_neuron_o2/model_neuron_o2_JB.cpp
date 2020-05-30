
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void JB_model_neuron_o2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JB[0+0*10] = x[0]*(-2.0/2.5E1)-5.0;
  JB[0+1*10] = -p[0]*p[1];
  JB[1+0*10] = 1.0;
  JB[1+1*10] = p[0];
  JB[2+0*10] = -x[2]*dwdx[1];
  JB[2+1*10] = -p[1];
  JB[2+2*10] = -w[1];
  JB[2+3*10] = -p[0]*p[1];
  JB[3+1*10] = 1.0;
  JB[3+2*10] = 1.0;
  JB[3+3*10] = p[0];
  JB[4+0*10] = -x[4]*dwdx[1];
  JB[4+1*10] = -p[0];
  JB[4+4*10] = -w[1];
  JB[4+5*10] = -p[0]*p[1];
  JB[5+4*10] = 1.0;
  JB[5+5*10] = p[0];
  JB[6+0*10] = -x[6]*dwdx[1];
  JB[6+6*10] = -w[1];
  JB[6+7*10] = -p[0]*p[1];
  JB[7+6*10] = 1.0;
  JB[7+7*10] = p[0];
  JB[8+0*10] = -x[8]*dwdx[1];
  JB[8+8*10] = -w[1];
  JB[8+9*10] = -p[0]*p[1];
  JB[9+8*10] = 1.0;
  JB[9+9*10] = p[0];
}

} // namespace model_model_neuron_o2

} // namespace amici

