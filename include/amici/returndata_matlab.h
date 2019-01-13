#ifndef RETURNDATA_MATLAB_H
#define RETURNDATA_MATLAB_H

#include "amici/rdata.h"

#include <vector>

#include <mex.h>

namespace amici {

mxArray *getReturnDataMatlabFromAmiciCall(ReturnData const *rdata);
mxArray *initMatlabReturnFields(ReturnData const *rdata);
mxArray *initMatlabDiagnosisFields(ReturnData const *rdata);

template <typename T>
void writeMatlabField0(mxArray *matlabStruct, const char *fieldName,
                       T fieldData);

template <typename T>
void writeMatlabField1(mxArray *matlabStruct, const char *fieldName,
                       std::vector<T> const &fieldData, const int dim0);

template <typename T>
void writeMatlabField2(mxArray *matlabStruct, const char *fieldName,
                       std::vector<T> const &fieldData, int dim0, int dim1,
                       std::vector<int> perm);

template <typename T>
void writeMatlabField3(mxArray *matlabStruct, const char *fieldName,
                       std::vector<T> const &fieldData, int dim0, int dim1,
                       int dim2, std::vector<int> perm);

template <typename T>
void writeMatlabField4(mxArray *matlabStruct, const char *fieldName,
                       std::vector<T> const &fieldData, int dim0, int dim1,
                       int dim2, int dim3, std::vector<int> perm);

double *initAndAttachArray(mxArray *matlabStruct, const char *fieldName,
                           std::vector<mwSize> dim);

void checkFieldNames(const char **fieldNames, const int fieldCount);

template <typename T>
std::vector<T> reorder(std::vector<T> const& input,
                       std::vector<int> const& order);

} // namespace amici

#endif // RETURNDATA_MATLAB_H
