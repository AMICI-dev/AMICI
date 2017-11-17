
#include <include/symbolic_functions.h>
#include <sundials/sundials_types.h> //realtype definition
#include <sundials/sundials_sparse.h> //SlsMat definition
#include <cmath> 

void JSparseB_model_nested_events(SlsMat JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 1;
  JSparseB->data[0] = p[4]-h[1]*p[3];
}

