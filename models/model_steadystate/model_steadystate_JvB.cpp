
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void JvB_model_steadystate(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx) {
  JvB[0] = vB[0]*(p[1]*x[1]+p[0]*dwdx[0]*2.0)+vB[1]*(p[1]*x[1]-p[0]*dwdx[0])-p[1]*x[1]*vB[2];
  JvB[1] = -vB[0]*(p[2]*2.0-p[1]*x[0])+vB[1]*(p[2]+p[1]*x[0])-p[1]*x[0]*vB[2];
  JvB[2] = vB[2]*(k[3]+dwdx[1])-vB[0]*dwdx[1]-vB[1]*dwdx[1];
}

