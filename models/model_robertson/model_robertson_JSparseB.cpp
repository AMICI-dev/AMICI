
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_robertson{

void JSparseB_model_robertson(SUNMatrixContent_Sparse JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexvals[1] = 1;
  JSparseB->indexvals[2] = 2;
  JSparseB->indexvals[3] = 0;
  JSparseB->indexvals[4] = 1;
  JSparseB->indexvals[5] = 2;
  JSparseB->indexvals[6] = 0;
  JSparseB->indexvals[7] = 1;
  JSparseB->indexvals[8] = 2;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 3;
  JSparseB->indexptrs[2] = 6;
  JSparseB->indexptrs[3] = 9;
  JSparseB->data[0] = cj+p[0];
  JSparseB->data[1] = -dwdx[0];
  JSparseB->data[2] = -dwdx[1];
  JSparseB->data[3] = -p[0];
  JSparseB->data[4] = cj+dwdx[0]+p[2]*x[1]*2.0;
  JSparseB->data[5] = dwdx[1];
  JSparseB->data[6] = -1.0;
  JSparseB->data[7] = -1.0;
  JSparseB->data[8] = -1.0;
}

} // namespace model_model_robertson

} // namespace amici

