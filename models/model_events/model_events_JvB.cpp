
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void JvB_model_events(realtype *JvB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *vB, const realtype *w, const realtype *dwdx) {
  JvB[0] = -p[1]*vB[1]*exp(t*(-1.0/1.0E1))+h[3]*p[0]*vB[0];
  JvB[1] = p[2]*vB[1];
  JvB[2] = vB[2];
}

