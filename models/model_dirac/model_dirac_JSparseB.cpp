
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_dirac{

void JSparseB_model_dirac(SUNMatrixContent_Sparse JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexvals[1] = 0;
  JSparseB->indexvals[2] = 1;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 1;
  JSparseB->indexptrs[2] = 3;
  JSparseB->data[0] = p[0];
  JSparseB->data[1] = -p[2];
  JSparseB->data[2] = p[3];
}

} // namespace model_model_dirac

} // namespace amici

