
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void JDiag_model_calvetti(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) {
  JDiag[0+0*6] = -cj-w[0]*(1.0E2/8.99E2)+dwdx[6];
  JDiag[1+0*6] = -cj-w[6]*1.202935161794779E-2+dwdx[19];
  JDiag[2+0*6] = -cj-w[31]*8.196729508204918E-9+dwdx[32];
  JDiag[3+0*6] = -dwdx[36]-w[23]*w[24]*w[25]*dwdx[34]-1.0;
  JDiag[4+0*6] = -dwdx[43]-1.0;
  JDiag[5+0*6] = -dwdx[52]-1.0;
}

} // namespace model_model_calvetti

} // namespace amici

