
#include <include/symbolic_functions.h>
#include <include/amici_defines.h> //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

void J_model_robertson(realtype *J, const realtype t, const realtype *x, const double *p, const double *k, const realtype *h, const realtype cj, const realtype *dx, const realtype *w, const realtype *dwdx) {
  J[0+0*3] = -cj-p[0];
  J[0+1*3] = dwdx[0];
  J[0+2*3] = dwdx[1];
  J[1+0*3] = p[0];
  J[1+1*3] = -cj-dwdx[0]-p[2]*x[1]*2.0;
  J[1+2*3] = -dwdx[1];
  J[2+0*3] = 1.0;
  J[2+1*3] = 1.0;
  J[2+2*3] = 1.0;
}

