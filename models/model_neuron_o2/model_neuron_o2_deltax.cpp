
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void deltax_model_neuron_o2(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old) {
              switch(ie) { 
              case 0: {
  deltax[0] = -p[2]-x[0];
  deltax[1] = p[3];
  deltax[2] = -(x[2]*(p[2]*5.0+p[3]+x[0]*5.0-(p[2]*p[2])*(1.0/2.5E1)+(x[0]*x[0])*(1.0/2.5E1)))/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  deltax[3] = -(x[2]*(p[0]*(p[3]+x[1]+p[1]*p[2])-p[0]*(x[1]-p[1]*x[0])))/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  deltax[4] = -(x[4]*(p[2]*5.0+p[3]+x[0]*5.0-(p[2]*p[2])*(1.0/2.5E1)+(x[0]*x[0])*(1.0/2.5E1)))/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  deltax[5] = -(x[4]*(p[0]*(p[3]+x[1]+p[1]*p[2])-p[0]*(x[1]-p[1]*x[0])))/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  deltax[6] = -(x[6]*(p[2]*5.0+p[3]+x[0]*5.0-(p[2]*p[2])*(1.0/2.5E1)+(x[0]*x[0])*(1.0/2.5E1)))/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2)-1.0;
  deltax[7] = -(x[6]*(p[0]*(p[3]+x[1]+p[1]*p[2])-p[0]*(x[1]-p[1]*x[0])))/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  deltax[8] = -(x[8]*(p[2]*5.0+p[3]+x[0]*5.0-(p[2]*p[2])*(1.0/2.5E1)+(x[0]*x[0])*(1.0/2.5E1)))/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2);
  deltax[9] = -(x[8]*(p[0]*(p[3]+x[1]+p[1]*p[2])-p[0]*(x[1]-p[1]*x[0])))/(k[1]+x[0]*5.0-x[1]+(x[0]*x[0])*(1.0/2.5E1)+1.4E2)+1.0;

              } break;

              } 
}

} // namespace model_model_neuron_o2

} // namespace amici

