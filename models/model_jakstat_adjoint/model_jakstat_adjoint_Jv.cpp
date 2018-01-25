
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void Jv_model_jakstat_adjoint(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx) {
  Jv[0] = -p[0]*v[0]*w[0]+(k[1]*p[3]*v[8])/k[0];
  Jv[1] = p[0]*v[0]*w[0]-p[1]*v[1]*dwdx[0]*2.0;
  Jv[2] = -p[2]*v[2]+p[1]*v[1]*dwdx[0];
  Jv[3] = -p[3]*v[3]+(k[0]*p[2]*v[2])/k[1];
  Jv[4] = p[3]*v[3]*2.0-p[3]*v[4];
  Jv[5] = p[3]*v[4]-p[3]*v[5];
  Jv[6] = p[3]*v[5]-p[3]*v[6];
  Jv[7] = p[3]*v[6]-p[3]*v[7];
  Jv[8] = p[3]*v[7]-p[3]*v[8];
}

