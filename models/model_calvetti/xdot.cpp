
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void xdot_model_calvetti(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *dx, const realtype *w) {
  xdot[0] = h[0]*(-1.0/3.1E1)+w[12]+w[22]-dx[0]-w[0]*x[0]*(1.0E2/8.99E2)+1.29E2/8.99E2;
  xdot[1] = w[30]-dx[1]-x[1]*w[6]*1.202935161794779E-2+1.816432094310117E-2;
  xdot[2] = w[37]-dx[2]-x[2]*w[31]*8.196729508204918E-9+8.196729508204918E-3;
  xdot[3] = h[0]*(1.0/3.1E1)-x[3]-w[12]-w[22]+w[0]*x[0]*(1.0E2/8.99E2)+w[23]*w[24]*w[25]*2.0-w[20]*w[23]*w[24]*w[25]-1.29E2/8.99E2;
  xdot[4] = x[3]-x[4]-w[30]+x[1]*w[6]*1.202935161794779E-2-1.816432094310117E-2;
  xdot[5] = x[4]-x[5]-w[37]+x[2]*w[31]*8.196729508204918E-9-8.196729508204918E-3;
}

} // namespace model_model_calvetti

} // namespace amici

