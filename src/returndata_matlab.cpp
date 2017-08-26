#include "returndata_matlab.h"

ReturnDataMatlab::ReturnDataMatlab(const UserData *udata, const Model *model)
    : ReturnData(udata, model) {
    mxsol = NULL;
    freeFieldsOnDestruction = false;
    initFields();
    copyFromUserData(udata);
}

void ReturnDataMatlab::initFields() {
    const int numFields = 35;
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
    mxArray *array = mxCreateDoubleMatrix(dim, 1, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if (status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}

void ReturnDataMatlab::initField2(double **fieldPointer, const char *fieldName,
                                  int dim1, int dim2) {
    mxArray *array = mxCreateDoubleMatrix(dim1, dim2, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if (status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}

void ReturnDataMatlab::initField3(double **fieldPointer, const char *fieldName,
                                  int dim1, int dim2, int dim3) {
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
    mwSize dims[] = {(mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3),
                     (mwSize)(dim4)};
    mxArray *array = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
    *fieldPointer = mxGetPr(array);
    mxSetField(mxsol, 0, fieldName, array);

    array = mxGetField(mxsol, 0, fieldName);
    if (status && array == NULL)
        *status = AMICI_ERROR_RDATA;
}
