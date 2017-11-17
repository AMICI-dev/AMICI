
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <cmath> 

void dfdx_model_robertson {
  dfdx[0+0*3] = -p[0];
  dfdx[0+1*3] = dwdx[0];
  dfdx[0+2*3] = dwdx[1];
  dfdx[1+0*3] = p[0];
  dfdx[1+1*3] = -dwdx[0]-p[2]*x[1]*2.0;
  dfdx[1+2*3] = -dwdx[1];
  dfdx[2+0*3] = 1.0;
  dfdx[2+1*3] = 1.0;
  dfdx[2+2*3] = 1.0;
}

