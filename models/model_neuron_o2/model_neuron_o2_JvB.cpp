
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void JvB_model_neuron_o2(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx) {
  JvB[0] = -vB[0]*(x[0]*(2.0/2.5E1)+5.0)-p[0]*p[1]*vB[1];
  JvB[1] = vB[0]+p[0]*vB[1];
  JvB[2] = -p[1]*vB[1]-w[1]*vB[2]-p[0]*p[1]*vB[3]-x[2]*vB[0]*dwdx[1];
  JvB[3] = vB[1]+vB[2]+p[0]*vB[3];
  JvB[4] = -p[0]*vB[1]-w[1]*vB[4]-p[0]*p[1]*vB[5]-x[4]*vB[0]*dwdx[1];
  JvB[5] = vB[4]+p[0]*vB[5];
  JvB[6] = -w[1]*vB[6]-p[0]*p[1]*vB[7]-x[6]*vB[0]*dwdx[1];
  JvB[7] = vB[6]+p[0]*vB[7];
  JvB[8] = -w[1]*vB[8]-p[0]*p[1]*vB[9]-x[8]*vB[0]*dwdx[1];
  JvB[9] = vB[8]+p[0]*vB[9];
}

