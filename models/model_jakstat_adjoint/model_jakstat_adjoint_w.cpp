
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void w_model_jakstat_adjoint(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h) {
  w[0] = amici::spline_pos(t,5,0.0,p[5],5.0,p[6],1.0E1,p[7],2.0E1,p[8],6.0E1,p[9],0.0,0.0);
  w[1] = x[1]*x[1];
}

