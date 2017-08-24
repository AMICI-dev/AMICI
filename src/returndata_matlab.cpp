#include "returndata_matlab.h"

ReturnDataMatlab::ReturnDataMatlab(const UserData *udata, const Model *model)
    : ReturnData(udata, model) {
        /**
         * @brief initialises the returnData struct, initialises the fields and copies
         * model dimensions from the udata struct
         * @param[in] udata pointer to the user data struct @type UserData
         * @param[in] model pointer to model specification object @type Model
         */
    mxsol = NULL;
    freeFieldsOnDestruction = false;
    initFields();
    copyFromUserData(udata);
}

void ReturnDataMatlab::initFields() {
    /**
     * @brief initialises sol object with the corresponding fields
     */
    const int numFields = 38;
    const char *field_names_sol[numFields] = {"status",
                                              "llh",
                                              "sllh",
                                              "s2llh",
                                              "chi2",
                                              "t",
                                              "x",
                                              "sx",
                                              "y",
                                              "sy",
                                              "sigmay",
                                              "ssigmay",
                                              "z",
                                              "sz",
                                              "sigmaz",
                                              "ssigmaz",
                                              "rz",
                                              "srz",
                                              "s2rz",
                                              "xdot",
                                              "J",
                                              "dydp",
                                              "dydx",
                                              "dxdotdp",
                                              "numsteps",
                                              "numrhsevals",
                                              "numerrtestfails",
                                              "numnonlinsolvconvfails",
                                              "order",
                                              "numstepsB",
                                              "numrhsevalsB",
                                              "numerrtestfailsB",
                                              "numnonlinsolvconvfailsB",
                                              "xss",
                                              "newton_status",
                                              "newton_numsteps",
                                              "newton_numlinsteps",
                                              "newton_time"};

    mxsol = mxCreateStructMatrix(1, 1, numFields, field_names_sol);

    ReturnData::initFields();
}

void ReturnDataMatlab::initField1(double **fieldPointer, const char *fieldName,
                                  int dim) {
    /**
     * @brief initialise vector and attach to the field
     * @param fieldPointer pointer of the field to which the vector will be attached
     * @param fieldName Name of the field to which the vector will be attached
     * @param dim number of elements in the vector
     */
    mxArray *array = mxCreateDoubleMatrix(dim, 1, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if (status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}

void ReturnDataMatlab::initField2(double **fieldPointer, const char *fieldName,
                                  int dim1, int dim2) {
    /**
     * @brief initialise matrix and attach to the field
     * @param fieldPointer pointer of the field to which the matrix will be attached
     * @param fieldName Name of the field to which the matrix will be attached
     * @param dim1 number of rows in the matrix
     * @param dim2 number of columns in the matrix
     */
    mxArray *array = mxCreateDoubleMatrix(dim1, dim2, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if (status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}

void ReturnDataMatlab::initField3(double **fieldPointer, const char *fieldName,
                                  int dim1, int dim2, int dim3) {
    /**
     * @brief initialise 3D tensor and attach to the field
     * @param fieldPointer pointer of the field to which the tensor will be attached
     * @param fieldName Name of the field to which the tensor will be attached
     * @param dim1 number of rows in the tensor
     * @param dim2 number of columns in the tensor
     * @param dim3 number of elements in the third dimension of the tensor
     */
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3)};
    mxArray *array = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if (status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}

void ReturnDataMatlab::initField4(double **fieldPointer, const char *fieldName,
                                  int dim1, int dim2, int dim3, int dim4) {
    /**
     * @brief initialise 4D tensor and attach to the field
     * @param fieldPointer pointer of the field to which the tensor will be attached
     * @param fieldName Name of the field to which the tensor will be attached
     * @param dim1 number of rows in the tensor
     * @param dim2 number of columns in the tensor
     * @param dim3 number of elements in the third dimension of the tensor
     * @param dim4 number of elements in the fourth dimension of the tensor
     */
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3),
                     (mwSize)(dim4)};
    mxArray *array = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if (status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}
