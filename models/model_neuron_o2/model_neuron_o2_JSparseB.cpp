
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron_o2{

void JSparseB_model_neuron_o2(SUNMatrixContent_Sparse JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexvals[1] = 1;
  JSparseB->indexvals[2] = 2;
  JSparseB->indexvals[3] = 4;
  JSparseB->indexvals[4] = 6;
  JSparseB->indexvals[5] = 8;
  JSparseB->indexvals[6] = 0;
  JSparseB->indexvals[7] = 1;
  JSparseB->indexvals[8] = 2;
  JSparseB->indexvals[9] = 3;
  JSparseB->indexvals[10] = 4;
  JSparseB->indexvals[11] = 2;
  JSparseB->indexvals[12] = 3;
  JSparseB->indexvals[13] = 2;
  JSparseB->indexvals[14] = 3;
  JSparseB->indexvals[15] = 4;
  JSparseB->indexvals[16] = 5;
  JSparseB->indexvals[17] = 4;
  JSparseB->indexvals[18] = 5;
  JSparseB->indexvals[19] = 6;
  JSparseB->indexvals[20] = 7;
  JSparseB->indexvals[21] = 6;
  JSparseB->indexvals[22] = 7;
  JSparseB->indexvals[23] = 8;
  JSparseB->indexvals[24] = 9;
  JSparseB->indexvals[25] = 8;
  JSparseB->indexvals[26] = 9;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 6;
  JSparseB->indexptrs[2] = 11;
  JSparseB->indexptrs[3] = 13;
  JSparseB->indexptrs[4] = 15;
  JSparseB->indexptrs[5] = 17;
  JSparseB->indexptrs[6] = 19;
  JSparseB->indexptrs[7] = 21;
  JSparseB->indexptrs[8] = 23;
  JSparseB->indexptrs[9] = 25;
  JSparseB->indexptrs[10] = 27;
  JSparseB->data[0] = x[0]*(-2.0/2.5E1)-5.0;
  JSparseB->data[1] = 1.0;
  JSparseB->data[2] = -x[2]*dwdx[1];
  JSparseB->data[3] = -x[4]*dwdx[1];
  JSparseB->data[4] = -x[6]*dwdx[1];
  JSparseB->data[5] = -x[8]*dwdx[1];
  JSparseB->data[6] = -p[0]*p[1];
  JSparseB->data[7] = p[0];
  JSparseB->data[8] = -p[1];
  JSparseB->data[9] = 1.0;
  JSparseB->data[10] = -p[0];
  JSparseB->data[11] = -w[1];
  JSparseB->data[12] = 1.0;
  JSparseB->data[13] = -p[0]*p[1];
  JSparseB->data[14] = p[0];
  JSparseB->data[15] = -w[1];
  JSparseB->data[16] = 1.0;
  JSparseB->data[17] = -p[0]*p[1];
  JSparseB->data[18] = p[0];
  JSparseB->data[19] = -w[1];
  JSparseB->data[20] = 1.0;
  JSparseB->data[21] = -p[0]*p[1];
  JSparseB->data[22] = p[0];
  JSparseB->data[23] = -w[1];
  JSparseB->data[24] = 1.0;
  JSparseB->data[25] = -p[0]*p[1];
  JSparseB->data[26] = p[0];
}

} // namespace model_model_neuron_o2

} // namespace amici

