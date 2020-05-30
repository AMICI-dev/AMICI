
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_nested_events{

void JSparseB_model_nested_events(SUNMatrixContent_Sparse JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 1;
  JSparseB->data[0] = p[4]+p[3]*(h[0]-1.0);
}

} // namespace model_model_nested_events

} // namespace amici

