#ifndef RETURNDATA_MATLAB_H
#define RETURNDATA_MATLAB_H

#include "include/rdata.h"
#include <mex.h>

class UserData;
class Model;

/**
 * @brief The ReturnDataMatlab class sets up ReturnData to be returned by the
 * MATLAB mex functions.
 * Memory is allocated using matlab functions.
 */
class ReturnDataMatlab : public ReturnData {

  public:
    ReturnDataMatlab(const UserData *udata, const Model *model);

    /** sol struct that is passed back to matlab */
    mxArray *matlabSolutionStruct = nullptr;

  protected:
    void initFields() override;

    virtual void initField1(double **fieldPointer, const char *fieldName,
                            int dim) override;

    virtual void initField2(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2) override;

    virtual void initField3(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2, int dim3) override;

    virtual void initField4(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2, int dim3, int dim4) override;
};

#endif // RETURNDATA_MATLAB_H
