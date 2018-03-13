
#include "amici/symbolic_functions.h"
#include "amici/amici_defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void Jv_model_steadystate(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *v, const realtype *w, const realtype *dwdx) {
  Jv[0] = v[2]*dwdx[1]+v[1]*(p[2]*2.0-p[1]*x[0])-v[0]*(p[1]*x[1]+p[0]*dwdx[0]*2.0);
  Jv[1] = v[2]*dwdx[1]-v[0]*(p[1]*x[1]-p[0]*dwdx[0])-v[1]*(p[2]+p[1]*x[0]);
  Jv[2] = -v[2]*(k[3]+dwdx[1])+p[1]*v[0]*x[1]+p[1]*v[1]*x[0];
}

