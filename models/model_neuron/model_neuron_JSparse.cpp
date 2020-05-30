
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
#include <sunmatrix/sunmatrix_sparse.h> //SUNMatrixContent_Sparse definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

namespace amici {

namespace model_model_neuron{

void JSparse_model_neuron(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx) {
  JSparse->indexvals[0] = 0;
  JSparse->indexvals[1] = 1;
  JSparse->indexvals[2] = 0;
  JSparse->indexvals[3] = 1;
  JSparse->indexptrs[0] = 0;
  JSparse->indexptrs[1] = 2;
  JSparse->indexptrs[2] = 4;
  JSparse->data[0] = x[0]*(2.0/2.5E1)+5.0;
  JSparse->data[1] = p[0]*p[1];
  JSparse->data[2] = -1.0;
  JSparse->data[3] = -p[0];
}

} // namespace model_model_neuron

} // namespace amici

