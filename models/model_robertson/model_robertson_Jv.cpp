
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void Jv_model_robertson(realtype *Jv, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *v, const realtype *w, const realtype *dwdx) {
  Jv[0] = -v[0]*(cj+p[0])+v[1]*dwdx[0]+v[2]*dwdx[1];
  Jv[1] = p[0]*v[0]-v[2]*dwdx[1]-v[1]*(cj+dwdx[0]+p[2]*x[1]*2.0);
  Jv[2] = v[0]+v[1]+v[2];
}

