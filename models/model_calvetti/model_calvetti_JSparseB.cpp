
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_calvetti{

void JSparseB_model_calvetti(SUNMatrixContent_Sparse JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype cj, const realtype *xB, const realtype *dx, const realtype *dxB, const realtype *w, const realtype *dwdx) {
  JSparseB->indexvals[0] = 0;
  JSparseB->indexvals[1] = 1;
  JSparseB->indexvals[2] = 2;
  JSparseB->indexvals[3] = 3;
  JSparseB->indexvals[4] = 4;
  JSparseB->indexvals[5] = 5;
  JSparseB->indexvals[6] = 1;
  JSparseB->indexvals[7] = 2;
  JSparseB->indexvals[8] = 4;
  JSparseB->indexvals[9] = 5;
  JSparseB->indexvals[10] = 2;
  JSparseB->indexvals[11] = 5;
  JSparseB->indexvals[12] = 0;
  JSparseB->indexvals[13] = 1;
  JSparseB->indexvals[14] = 2;
  JSparseB->indexvals[15] = 3;
  JSparseB->indexvals[16] = 4;
  JSparseB->indexvals[17] = 5;
  JSparseB->indexvals[18] = 1;
  JSparseB->indexvals[19] = 2;
  JSparseB->indexvals[20] = 3;
  JSparseB->indexvals[21] = 4;
  JSparseB->indexvals[22] = 5;
  JSparseB->indexvals[23] = 2;
  JSparseB->indexvals[24] = 4;
  JSparseB->indexvals[25] = 5;
  JSparseB->indexptrs[0] = 0;
  JSparseB->indexptrs[1] = 6;
  JSparseB->indexptrs[2] = 10;
  JSparseB->indexptrs[3] = 12;
  JSparseB->indexptrs[4] = 18;
  JSparseB->indexptrs[5] = 23;
  JSparseB->indexptrs[6] = 26;
  JSparseB->data[0] = cj+w[0]*(1.0E2/8.99E2)-dwdx[6];
  JSparseB->data[1] = -dwdx[16];
  JSparseB->data[2] = -dwdx[28];
  JSparseB->data[3] = -dwdx[36];
  JSparseB->data[4] = -dwdx[40];
  JSparseB->data[5] = -dwdx[47];
  JSparseB->data[6] = cj+w[6]*1.202935161794779E-2-dwdx[19];
  JSparseB->data[7] = -dwdx[31];
  JSparseB->data[8] = -dwdx[43];
  JSparseB->data[9] = -dwdx[50];
  JSparseB->data[10] = cj+w[31]*8.196729508204918E-9-dwdx[32];
  JSparseB->data[11] = -dwdx[52];
  JSparseB->data[12] = w[0]*(-1.0E2/8.99E2)+dwdx[6]-dwdx[7]*(w[23]*w[24]*2.0-w[20]*w[23]*w[24])+w[23]*w[24]*w[25]*dwdx[4];
  JSparseB->data[13] = dwdx[16]+w[23]*w[24]*w[25]*dwdx[14];
  JSparseB->data[14] = dwdx[28]+w[23]*w[24]*w[25]*dwdx[26];
  JSparseB->data[15] = dwdx[36]+w[23]*w[24]*w[25]*dwdx[34]+1.0;
  JSparseB->data[16] = dwdx[40]+w[23]*w[24]*w[25]*dwdx[38];
  JSparseB->data[17] = dwdx[47]+w[23]*w[24]*w[25]*dwdx[45];
  JSparseB->data[18] = w[6]*(-1.202935161794779E-2)+dwdx[19];
  JSparseB->data[19] = dwdx[31];
  JSparseB->data[20] = -1.0;
  JSparseB->data[21] = dwdx[43]+1.0;
  JSparseB->data[22] = dwdx[50];
  JSparseB->data[23] = w[31]*(-8.196729508204918E-9)+dwdx[32];
  JSparseB->data[24] = -1.0;
  JSparseB->data[25] = dwdx[52]+1.0;
}

} // namespace model_model_calvetti

} // namespace amici

