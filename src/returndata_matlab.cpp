#include "returndata_matlab.h"

#include <rdata.h>

namespace amici {

mxArray *getReturnDataMatlabFromAmiciCall(ReturnData const *rdata) {
    /**
      * @brief generates matlab mxArray from a ReturnData object
      * @param rdata ReturnDataObject
      * @return rdatamatlab ReturnDataObject stored as matlab compatible data
      */
    
    mxArray *matlabSolutionStruct = initMatlabReturnFields(rdata);
    return matlabSolutionStruct;
}

mxArray *initMatlabReturnFields(ReturnData const *rdata) {
    /**
     * @brief initialises sol object with the corresponding fields
     * @param rdata ReturnDataObject
     */
    const int numFields = 22;
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
                                              "x0",
                                              "sx0",
                                              "diagnosis"};

    mxArray *matlabSolutionStruct =
        mxCreateStructMatrix(1, 1, numFields, field_names_sol);
    
    writeMatlabField0(matlabSolutionStruct, "status", rdata->status);

    writeMatlabField1(matlabSolutionStruct, "t", rdata->ts, rdata->nt);
    writeMatlabField0(matlabSolutionStruct, "llh", rdata->llh);
    writeMatlabField0(matlabSolutionStruct, "chi2", rdata->chi2);

    if ((rdata->nz > 0) & (rdata->ne > 0)) {
        writeMatlabField2(matlabSolutionStruct, "z", rdata->z, rdata->nmaxevent, rdata->nz);
        writeMatlabField2(matlabSolutionStruct, "rz", rdata->rz, rdata->nmaxevent, rdata->nz);
        writeMatlabField2(matlabSolutionStruct, "sigmaz", rdata->sigmaz, rdata->nmaxevent, rdata->nz);
    }
    if (rdata->nx > 0) {
        writeMatlabField2(matlabSolutionStruct, "x0",  rdata->x0, 1, rdata->nx);
        writeMatlabField2(matlabSolutionStruct, "sx0",  rdata->sx0, rdata->nx, rdata->nplist);
    }
    if (rdata->ny > 0) {
        writeMatlabField2(matlabSolutionStruct, "y", rdata->y, rdata->nt, rdata->ny);
        writeMatlabField2(matlabSolutionStruct, "sigmay", rdata->sigmay, rdata->nt, rdata->ny);
    }
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        writeMatlabField2(matlabSolutionStruct, "sllh", rdata->sllh, rdata->nplist, 1);

        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            writeMatlabField3(matlabSolutionStruct, "sx", rdata->sx, rdata->nt, rdata->nx, rdata->nplist);
            if (rdata->ny > 0) {
                writeMatlabField3(matlabSolutionStruct, "sy", rdata->sy, rdata->nt, rdata->ny, rdata->nplist);
            }
            if ((rdata->nz > 0) & (rdata->ne > 0)) {
                writeMatlabField3(matlabSolutionStruct, "srz", rdata->srz, rdata->nmaxevent, rdata->nz, rdata->nplist);
                if (rdata->sensi >= AMICI_SENSI_ORDER_SECOND) {
                    writeMatlabField4(matlabSolutionStruct, "s2rz", rdata->s2rz, rdata->nmaxevent, rdata->nztrue, rdata->nplist,
                               rdata->nplist);
                }
                writeMatlabField3(matlabSolutionStruct, "sz", rdata->sz, rdata->nmaxevent, rdata->nz, rdata->nplist);
            }
        }

        if (rdata->ny > 0) {
            writeMatlabField3(matlabSolutionStruct, "ssigmay", rdata->ssigmay, rdata->nt, rdata->ny, rdata->nplist);
        }
        if ((rdata->nz > 0) & (rdata->ne > 0)) {
            writeMatlabField3(matlabSolutionStruct, "ssigmaz", rdata->ssigmaz, rdata->nmaxevent, rdata->nz, rdata->nplist);
        }

        if (rdata->sensi >= AMICI_SENSI_ORDER_SECOND) {
            writeMatlabField2(matlabSolutionStruct, "s2llh", rdata->s2llh, rdata->nJ - 1, rdata->nplist);
        }
    }

    mxArray *diagnosis = initMatlabDiagnosisFields(rdata);
    mxSetField(matlabSolutionStruct, 0, "diagnosis", diagnosis);
    
    return(matlabSolutionStruct);
}

mxArray *initMatlabDiagnosisFields(ReturnData const *rdata) {
    /**
     * @brief initialises diagnosis object with the corresponding fields
     * @param rdata ReturnDataObject
     */
    const int numFields = 15;
    const char *field_names_sol[numFields] = {"xdot",
                                              "J",
                                              "numsteps",
                                              "numrhsevals",
                                              "numerrtestfails",
                                              "numnonlinsolvconvfails",
                                              "order",
                                              "numstepsB",
                                              "numrhsevalsB",
                                              "numerrtestfailsB",
                                              "numnonlinsolvconvfailsB"
                                              "newton_status",
                                              "newton_numsteps",
                                              "newton_numlinsteps",
                                              "newton_time"};

    mxArray *matlabDiagnosisStruct =
        mxCreateStructMatrix(1, 1, numFields, field_names_sol);
    
    writeMatlabField1(matlabDiagnosisStruct, "numsteps", rdata->numsteps, rdata->nt);
    writeMatlabField1(matlabDiagnosisStruct, "numrhsevals", rdata->numrhsevals, rdata->nt);
    writeMatlabField1(matlabDiagnosisStruct, "numerrtestfails", rdata->numrhsevals, rdata->nt);
    writeMatlabField1(matlabDiagnosisStruct, "numnonlinsolvconvfails", rdata->numnonlinsolvconvfails, rdata->nt);
    writeMatlabField1(matlabDiagnosisStruct, "order", rdata->order, rdata->nt);

    if (rdata->nx > 0) {
        writeMatlabField1(matlabDiagnosisStruct, "xdot", rdata->xdot, rdata->nx);
        writeMatlabField2(matlabDiagnosisStruct, "J", rdata->J, rdata->nx, rdata->nx);
        writeMatlabField0(matlabDiagnosisStruct, "newton_status", rdata->newton_status);
        writeMatlabField1(matlabDiagnosisStruct, "newton_numsteps", rdata->newton_numsteps, 2);
        writeMatlabField2(matlabDiagnosisStruct, "newton_numlinsteps", rdata->newton_numlinsteps, rdata->newton_maxsteps, 2);
        writeMatlabField0(matlabDiagnosisStruct, "newton_time", rdata->newton_time);
    }
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        if (rdata->sensi_meth == AMICI_SENSI_ASA) {
            writeMatlabField1(matlabDiagnosisStruct, "numstepsB", rdata->numstepsB, rdata->nt);
            writeMatlabField1(matlabDiagnosisStruct, "numrhsevalsB", rdata->numrhsevalsB, rdata->nt);
            writeMatlabField1(matlabDiagnosisStruct, "numerrtestfailsB", rdata->numrhsevals, rdata->nt);
            writeMatlabField1(matlabDiagnosisStruct, "numnonlinsolvconvfailsB", rdata->numnonlinsolvconvfailsB, rdata->nt);
        }
    }

    return(matlabDiagnosisStruct);
}


template<typename T>
void writeMatlabField0(mxArray *matlabSolutionStruct, const char *fieldName,
                       T fielddata) {
    /**
     * @brief initialise vector and attach to the field
     * @param fieldPointer pointer of the field to which the vector will be
     * attached
     * @param fieldName Name of the field to which the vector will be attached
     */
    
    mxArray *matlab_array = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxSetField(matlabSolutionStruct, 0, fieldName, matlab_array);
    double *c_array = mxGetPr(matlab_array);
    c_array[0] = static_cast<double>(fielddata);
}

template<typename T>
void writeMatlabField1(mxArray *matlabSolutionStruct, const char *fieldName,
                       const std::vector<T> fieldData, int dim0) {
    /**
     * @brief initialise vector and attach to the field
     * @param fieldPointer pointer of the field to which the vector will be
     * attached
     * @param fieldName Name of the field to which the vector will be attached
     * @param dim0 number of elements in the vector
     */
    if(fieldData.size() != dim0) {
        throw AmiException("Dimension mismatch when writing rdata->%s to matlab results",fieldName);
    }
    
    mxArray *matlab_array = mxCreateDoubleMatrix(dim0, 1, mxREAL);
    mxSetField(matlabSolutionStruct, 0, fieldName, matlab_array);
    double *c_array = mxGetPr(matlab_array);
    for(int i = 0; i < dim0; i++)
        c_array[i] = static_cast<double>(fieldData[i]);
}

template<typename T>
void writeMatlabField2(mxArray *matlabStruct, const char *fieldName,
                      const std::vector<T> fieldData, int dim0, int dim1) {
    /**
     * @brief initialise matrix, attach to the field and write data
     * @param matlabStruct Pointer to the matlab structure
     * @param fieldName Name of the field to which the tensor will be attached
     * @param fieldData Data wich will be stored in the field
     * @param dim0 number of rows in the tensor
     * @param dim1 number of columns in the tensor
     */
    if(fieldData.size() != dim0*dim1) {
        throw AmiException("Dimension mismatch when writing rdata->%s to matlab results",fieldName);
    }
    std::vector<mwSize> dim = {(mwSize)(dim0), (mwSize)(dim1)};
    
    double *array = initAndAttachArray(matlabStruct, fieldName, dim);
    
    /* transform rowmajor (c++) to colmajor (matlab) */
    for (int i = 0; i < dim0; i++) {
        for (int j = 0; j < dim1; j++) {
            array[i + j*dim0] =
                static_cast<double>(fieldData[i*dim1 + j]);
        }
    }
}

template<typename T>
void writeMatlabField3(mxArray *matlabStruct, const char *fieldName,
                      const std::vector<T> fieldData, int dim0, int dim1,
                      int dim2) {
    /**
     * @brief initialise 3D tensor, attach to the field and write data
     * @param matlabStruct Pointer to the matlab structure
     * @param fieldName Name of the field to which the tensor will be attached
     * @param fieldData Data wich will be stored in the field
     * @param dim0 number of rows in the tensor
     * @param dim1 number of columns in the tensor
     * @param dim2 number of elements in the third dimension of the tensor
     */
    if(fieldData.size() != dim0*dim1*dim2) {
        throw AmiException("Dimension mismatch when writing rdata->%s to matlab results",fieldName);
    }
    
    std::vector<mwSize> dim = {(mwSize)(dim0), (mwSize)(dim1), (mwSize)(dim2)};
    
    double *array = initAndAttachArray(matlabStruct, fieldName, dim);
    
    /* transform rowmajor (c++) to colmajor (matlab) */
    for (int i = 0; i < dim0; i++) {
        for (int j = 0; j < dim1; j++) {
            for (int k = 0; k < dim2; k++) {
                array[i + (j + k*dim1)*dim0] =
                static_cast<double>(fieldData[(i*dim1 + j)*dim2 + k]);
            }
        }
    }
}

template<typename T>
void writeMatlabField4(mxArray *matlabStruct, const char *fieldName,
                      const std::vector<T> fieldData, int dim0, int dim1,
                      int dim2, int dim3) {
    /**
     * @brief initialise 4D tensor, attach to the field and write data
     * @param matlabStruct Pointer to the matlab structure
     * @param fieldName Name of the field to which the tensor will be attached
     * @param fieldData Data wich will be stored in the field
     * @param dim0 number of rows in the tensor
     * @param dim1 number of columns in the tensor
     * @param dim2 number of elements in the third dimension of the tensor
     * @param dim3 number of elements in the fourth dimension of the tensor
     */
    if(fieldData.size() != dim0*dim1*dim2*dim3) {
        throw AmiException("Dimension mismatch when writing rdata->%s to matlab results",fieldName);
    }
    
    std::vector<mwSize> dim = {(mwSize)(dim0), (mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3)};
    
    double *array = initAndAttachArray(matlabStruct, fieldName, dim);
    
    /* transform rowmajor (c++) to colmajor (matlab) */
    for (int i = 0; i < dim0; i++) {
        for (int j = 0; j < dim1; j++) {
            for (int k = 0; k < dim2; k++) {
                for (int l = 0; l < dim3; l++) {
                    array[i + (j + (k + l*dim2)*dim1)*dim0] =
                    static_cast<double>(fieldData[((i*dim1 + j)*dim2 + k)*dim3 + l]);
                }
            }
        }
    }
}

double *initAndAttachArray(mxArray *matlabStruct, const char *fieldName, std::vector<mwSize> dim) {
    mxArray *array = mxCreateNumericArray(dim.size(), dim.data(), mxDOUBLE_CLASS, mxREAL);
    mxSetField(matlabStruct, 0, fieldName, array);
    return(mxGetPr(array));
}


} // namespace amici
