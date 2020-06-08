
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void JB_model_calvetti(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) {
  JB[0+0*6] = cj+w[0]*(1.0E2/8.99E2)-dwdx[6];
  JB[0+3*6] = w[0]*(-1.0E2/8.99E2)+dwdx[6]-dwdx[7]*(w[23]*w[24]*2.0-w[20]*w[23]*w[24])+w[23]*w[24]*w[25]*dwdx[4];
  JB[1+0*6] = -dwdx[16];
  JB[1+1*6] = cj+w[6]*1.202935161794779E-2-dwdx[19];
  JB[1+3*6] = dwdx[16]+w[23]*w[24]*w[25]*dwdx[14];
  JB[1+4*6] = w[6]*(-1.202935161794779E-2)+dwdx[19];
  JB[2+0*6] = -dwdx[28];
  JB[2+1*6] = -dwdx[31];
  JB[2+2*6] = cj+w[31]*8.196729508204918E-9-dwdx[32];
  JB[2+3*6] = dwdx[28]+w[23]*w[24]*w[25]*dwdx[26];
  JB[2+4*6] = dwdx[31];
  JB[2+5*6] = w[31]*(-8.196729508204918E-9)+dwdx[32];
  JB[3+0*6] = -dwdx[36];
  JB[3+3*6] = dwdx[36]+w[23]*w[24]*w[25]*dwdx[34]+1.0;
  JB[3+4*6] = -1.0;
  JB[4+0*6] = -dwdx[40];
  JB[4+1*6] = -dwdx[43];
  JB[4+3*6] = dwdx[40]+w[23]*w[24]*w[25]*dwdx[38];
  JB[4+4*6] = dwdx[43]+1.0;
  JB[4+5*6] = -1.0;
  JB[5+0*6] = -dwdx[47];
  JB[5+1*6] = -dwdx[50];
  JB[5+2*6] = -dwdx[52];
  JB[5+3*6] = dwdx[47]+w[23]*w[24]*w[25]*dwdx[45];
  JB[5+4*6] = dwdx[50];
  JB[5+5*6] = dwdx[52]+1.0;
}

} // namespace model_model_calvetti

} // namespace amici

