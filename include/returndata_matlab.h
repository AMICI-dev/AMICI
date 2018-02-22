#ifndef RETURNDATA_MATLAB_H
#define RETURNDATA_MATLAB_H

#include "include/rdata.h"
#include <mex.h>

namespace amici {

void initFields();

mxArray *getReturnDataMatlabFromAmiciCall(ReturnData const *rdata);
mxArray *initMatlabReturnFields(ReturnData const *rdata);
mxArray *initMatlabDiagnosisFields(ReturnData const *rdata);

template<typename T>
void writeMatlabField0(mxArray *matlabSolutionStruct, const char *fieldName, T fielddata);

template <typename T>
void writeMatlabField1(mxArray *matlabstruct, const char *fieldName, std::vector<T> fieldData, const int dim0);

template <typename T>
void writeMatlabField2(mxArray *matlabstruct, const char *fieldName, std::vector<T> fieldData, int dim0, int dim1, std::vector<int> perm);

template <typename T>
void writeMatlabField3(mxArray *matlabstruct, const char *fieldName, std::vector<T> fieldData, int dim0, int dim1, int dim2, std::vector<int> perm);

template <typename T>
void writeMatlabField4(mxArray *matlabstruct, const char *fieldName, std::vector<T> fieldData, int dim0, int dim1, int dim2, int dim3, std::vector<int> perm);

double *initAndAttachArray(mxArray *matlabStruct, const char *fieldName, std::vector<mwSize> dim);

void checkFieldNames(const char **fieldNames,const int fieldCount);

template<typename T>
std::vector<T> reorder(const std::vector<T> input, const std::vector<int> order);

} // namespace amici

#endif // RETURNDATA_MATLAB_H
