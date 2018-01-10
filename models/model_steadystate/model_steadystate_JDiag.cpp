
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void JDiag_model_steadystate(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JDiag[0+0*3] = -p[1]*x[1]-p[0]*dwdx[0]*2.0;
  JDiag[1+0*3] = -p[2]-p[1]*x[0];
  JDiag[2+0*3] = -k[3]-dwdx[1];
}

