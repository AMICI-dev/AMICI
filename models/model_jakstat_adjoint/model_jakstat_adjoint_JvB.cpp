
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void JvB_model_jakstat_adjoint(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx) {
  JvB[0] = p[0]*w[0]*vB[0]-p[0]*w[0]*vB[1];
  JvB[1] = p[1]*vB[1]*dwdx[0]*2.0-p[1]*vB[2]*dwdx[0];
  JvB[2] = p[2]*vB[2]-(k[0]*p[2]*vB[3])/k[1];
  JvB[3] = p[3]*vB[3]-p[3]*vB[4]*2.0;
  JvB[4] = p[3]*vB[4]-p[3]*vB[5];
  JvB[5] = p[3]*vB[5]-p[3]*vB[6];
  JvB[6] = p[3]*vB[6]-p[3]*vB[7];
  JvB[7] = p[3]*vB[7]-p[3]*vB[8];
  JvB[8] = p[3]*vB[8]-(k[1]*p[3]*vB[0])/k[0];
}

