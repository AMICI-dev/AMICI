
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron{

void JSparseB_model_neuron(SUNMatrixContent_Sparse JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexvals[1] = 1;
  JSparseB->indexvals[2] = 0;
  JSparseB->indexvals[3] = 1;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 2;
  JSparseB->indexptrs[2] = 4;
  JSparseB->data[0] = x[0]*(-2.0/2.5E1)-5.0;
  JSparseB->data[1] = 1.0;
  JSparseB->data[2] = -p[0]*p[1];
  JSparseB->data[3] = p[0];
}

} // namespace model_model_neuron

} // namespace amici

