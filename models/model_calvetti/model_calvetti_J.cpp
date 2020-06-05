
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void J_model_calvetti(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) {
  J[0+0*6] = -cj-w[0]*(1.0E2/8.99E2)+dwdx[6];
  J[0+1*6] = dwdx[16];
  J[0+2*6] = dwdx[28];
  J[0+3*6] = dwdx[36];
  J[0+4*6] = dwdx[40];
  J[0+5*6] = dwdx[47];
  J[1+1*6] = -cj-w[6]*1.202935161794779E-2+dwdx[19];
  J[1+2*6] = dwdx[31];
  J[1+4*6] = dwdx[43];
  J[1+5*6] = dwdx[50];
  J[2+2*6] = -cj-w[31]*8.196729508204918E-9+dwdx[32];
  J[2+5*6] = dwdx[52];
  J[3+0*6] = w[0]*(1.0E2/8.99E2)-dwdx[6]+dwdx[7]*(w[23]*w[24]*2.0-w[20]*w[23]*w[24])-w[23]*w[24]*w[25]*dwdx[4];
  J[3+1*6] = -dwdx[16]-w[23]*w[24]*w[25]*dwdx[14];
  J[3+2*6] = -dwdx[28]-w[23]*w[24]*w[25]*dwdx[26];
  J[3+3*6] = -dwdx[36]-w[23]*w[24]*w[25]*dwdx[34]-1.0;
  J[3+4*6] = -dwdx[40]-w[23]*w[24]*w[25]*dwdx[38];
  J[3+5*6] = -dwdx[47]-w[23]*w[24]*w[25]*dwdx[45];
  J[4+1*6] = w[6]*1.202935161794779E-2-dwdx[19];
  J[4+2*6] = -dwdx[31];
  J[4+3*6] = 1.0;
  J[4+4*6] = -dwdx[43]-1.0;
  J[4+5*6] = -dwdx[50];
  J[5+2*6] = w[31]*8.196729508204918E-9-dwdx[32];
  J[5+4*6] = 1.0;
  J[5+5*6] = -dwdx[52]-1.0;
}

} // namespace model_model_calvetti

} // namespace amici

