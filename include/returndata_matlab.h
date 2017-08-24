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
    ~ReturnDataMatlab() {}

    mxArray *mxsol;

  protected:
    void initFields();

    virtual void initField1(double **fieldPointer, const char *fieldName,
                            int dim);

    /**
     * @ brief initialise matrix and attach to the field
     * @ param FIELD name of the field to which the matrix will be attached
     * @ param D1 number of rows in the matrix
     * @ param D2 number of columns in the matrix
     */

    virtual void initField2(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2);

    /**
     * @ brief initialise 3D tensor and attach to the field
     * @ param FIELD name of the field to which the tensor will be attached
     * @ param D1 number of rows in the tensor
     * @ param D2 number of columns in the tensor
     * @ param D3 number of elements in the third dimension of the tensor
     */

    virtual void initField3(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2, int dim3);

    /**
     * @ brief initialise 4D tensor and attach to the field
     * @ param FIELD name of the field to which the tensor will be attached
     * @ param D1 number of rows in the tensor
     * @ param D2 number of columns in the tensor
     * @ param D3 number of elements in the third dimension of the tensor
     * @ param D4 number of elements in the fourth dimension of the tensor
     */

    virtual void initField4(double **fieldPointer, const char *fieldName,
                            int dim1, int dim2, int dim3, int dim4);
};

#endif // RETURNDATA_MATLAB_H
