
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sundials/sundials_sparse.h> //SlsMat definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void JSparse_model_nested_events(SlsMat JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JSparse->indexvals[0] = 0;
  JSparse->indexptrs[0] = 0;
  JSparse->indexptrs[1] = 1;
  JSparse->data[0] = -p[4]-p[3]*(h[0]-1.0);
}

