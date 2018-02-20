#ifndef RETURNDATA_MATLAB_H
#define RETURNDATA_MATLAB_H

#include "include/rdata.h"
#include <mex.h>

namespace amici {

void initFields();

template <typename T>
void writeMatlabField1(mxArray *matlabstruct, const char *fieldName, std::vector<T> fieldData, const int dim);

template <typename T>
void writeMatlabField2(mxArray *matlabstruct, const char *fieldName, std::vector<T> fieldData, int dim1, int dim2);

template <typename T>
void writeMatlabField3(mxArray *matlabstruct, const char *fieldName, std::vector<T> fieldData, int dim1, int dim2, int dim3);

template <typename T>
void writeMatlabField4(mxArray *matlabstruct, const char *fieldName, std::vector<T> fieldData, int dim1, int dim2, int dim3, int dim4);
};

} // namespace amici

#endif // RETURNDATA_MATLAB_H
