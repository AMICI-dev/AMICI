#include "amici/returndata_matlab.h"
#include "amici/exception.h"
#include "amici/defines.h"

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
    
    checkFieldNames(field_names_sol,numFields);

    mxArray *matlabSolutionStruct =
        mxCreateStructMatrix(1, 1, numFields, field_names_sol);
    
    std::vector<int> perm1 = {0, 1};
    std::vector<int> perm2 = {0, 2, 1};
    std::vector<int> perm3 = {0, 2, 3, 1};
    
    writeMatlabField0(matlabSolutionStruct, "status", rdata->status);

    writeMatlabField1(matlabSolutionStruct, "t", rdata->ts, rdata->nt);
    writeMatlabField0(matlabSolutionStruct, "llh", rdata->llh);
    writeMatlabField0(matlabSolutionStruct, "chi2", rdata->chi2);

    if ((rdata->nz > 0) & (rdata->ne > 0)) {
        writeMatlabField2(matlabSolutionStruct, "z", rdata->z, rdata->nmaxevent, rdata->nz, perm1);
        writeMatlabField2(matlabSolutionStruct, "rz", rdata->rz, rdata->nmaxevent, rdata->nz, perm1);
        writeMatlabField2(matlabSolutionStruct, "sigmaz", rdata->sigmaz, rdata->nmaxevent, rdata->nz, perm1);
    }
    if (rdata->nx > 0) {
        writeMatlabField2(matlabSolutionStruct, "x", rdata->x, rdata->nt, rdata->nx, perm1);
        writeMatlabField2(matlabSolutionStruct, "x0",  rdata->x0, 1, rdata->nx, perm1);
    }
    if (rdata->ny > 0) {
        writeMatlabField2(matlabSolutionStruct, "y", rdata->y, rdata->nt, rdata->ny, perm1);
        writeMatlabField2(matlabSolutionStruct, "sigmay", rdata->sigmay, rdata->nt, rdata->ny, perm1);
    }
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        writeMatlabField1(matlabSolutionStruct, "sllh", rdata->sllh, rdata->nplist);
        writeMatlabField2(matlabSolutionStruct, "sx0",  rdata->sx0, rdata->nx, rdata->nplist, perm1);

        if (rdata->sensi_meth == AMICI_SENSI_FSA) {
            writeMatlabField3(matlabSolutionStruct, "sx", rdata->sx, rdata->nt, rdata->nplist, rdata->nx, perm2);
            if (rdata->ny > 0) {
                writeMatlabField3(matlabSolutionStruct, "sy", rdata->sy, rdata->nt, rdata->nplist, rdata->ny, perm2);
            }
            if ((rdata->nz > 0) & (rdata->ne > 0)) {
                writeMatlabField3(matlabSolutionStruct, "srz", rdata->srz, rdata->nmaxevent, rdata->nplist, rdata->nz, perm2);
                if (rdata->sensi >= AMICI_SENSI_ORDER_SECOND) {
                    writeMatlabField4(matlabSolutionStruct, "s2rz", rdata->s2rz, rdata->nmaxevent, rdata->nplist, rdata->nztrue,
                               rdata->nplist, perm3);
                }
                writeMatlabField3(matlabSolutionStruct, "sz", rdata->sz, rdata->nmaxevent, rdata->nplist, rdata->nz, perm2);
            }
        }

        if (rdata->ny > 0) {
            writeMatlabField3(matlabSolutionStruct, "ssigmay", rdata->ssigmay, rdata->nt, rdata->nplist, rdata->ny, perm2);
        }
        if ((rdata->nz > 0) & (rdata->ne > 0)) {
            writeMatlabField3(matlabSolutionStruct, "ssigmaz", rdata->ssigmaz, rdata->nmaxevent, rdata->nplist, rdata->nz, perm2);
        }

        if (rdata->sensi >= AMICI_SENSI_ORDER_SECOND) {
            writeMatlabField2(matlabSolutionStruct, "s2llh", rdata->s2llh, rdata->nplist, rdata->nJ - 1, perm1);
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
                                              "numnonlinsolvconvfailsB",
                                              "newton_status",
                                              "newton_numsteps",
                                              "newton_numlinsteps",
                                              "newton_time"};
    
    checkFieldNames(field_names_sol,numFields);
    
    mxArray *matlabDiagnosisStruct =
        mxCreateStructMatrix(1, 1, numFields, field_names_sol);
    
    std::vector<int> perm1 = {0, 1};
    
    writeMatlabField1(matlabDiagnosisStruct, "numsteps", rdata->numsteps, rdata->nt);
    writeMatlabField1(matlabDiagnosisStruct, "numrhsevals", rdata->numrhsevals, rdata->nt);
    writeMatlabField1(matlabDiagnosisStruct, "numerrtestfails", rdata->numerrtestfails, rdata->nt);
    writeMatlabField1(matlabDiagnosisStruct, "numnonlinsolvconvfails", rdata->numnonlinsolvconvfails, rdata->nt);
    writeMatlabField1(matlabDiagnosisStruct, "order", rdata->order, rdata->nt);

    if (rdata->nx > 0) {
        writeMatlabField1(matlabDiagnosisStruct, "xdot", rdata->xdot, rdata->nx);
        writeMatlabField2(matlabDiagnosisStruct, "J", rdata->J, rdata->nx, rdata->nx, perm1);
        writeMatlabField0(matlabDiagnosisStruct, "newton_status", rdata->newton_status);
        writeMatlabField1(matlabDiagnosisStruct, "newton_numsteps", rdata->newton_numsteps, 2);
        writeMatlabField2(matlabDiagnosisStruct, "newton_numlinsteps", rdata->newton_numlinsteps, rdata->newton_maxsteps, 2, perm1);
        writeMatlabField0(matlabDiagnosisStruct, "newton_time", rdata->newton_time);
    }
    if (rdata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        if (rdata->sensi_meth == AMICI_SENSI_ASA) {
            writeMatlabField1(matlabDiagnosisStruct, "numstepsB", rdata->numstepsB, rdata->nt);
            writeMatlabField1(matlabDiagnosisStruct, "numrhsevalsB", rdata->numrhsevalsB, rdata->nt);
            writeMatlabField1(matlabDiagnosisStruct, "numerrtestfailsB", rdata->numerrtestfailsB, rdata->nt);
            writeMatlabField1(matlabDiagnosisStruct, "numnonlinsolvconvfailsB", rdata->numnonlinsolvconvfailsB, rdata->nt);
        }
    }

    return(matlabDiagnosisStruct);
}


template<typename T>
void writeMatlabField0(mxArray *matlabStruct, const char *fieldName,
                       T fielddata) {
    /**
     * @brief initialise vector and attach to the field
     * @param fieldPointer pointer of the field to which the vector will be
     * attached
     * @param fieldName Name of the field to which the vector will be attached
     */
    
    std::vector<mwSize> dim = {(mwSize)(1), (mwSize)(1)};
    
    double *array = initAndAttachArray(matlabStruct, fieldName, dim);
    
    array[0] = static_cast<double>(fielddata);
}

template<typename T>
void writeMatlabField1(mxArray *matlabStruct, const char *fieldName,
                       const std::vector<T> fieldData, int dim0) {
    /**
     * @brief initialise vector and attach to the field
     * @param fieldPointer pointer of the field to which the vector will be
     * attached
     * @param fieldName Name of the field to which the vector will be attached
     * @param dim0 number of elements in the vector
     */
    if(fieldData.size() != dim0)
        throw AmiException("Dimension mismatch when writing rdata->%s to matlab results",fieldName);
    
    std::vector<mwSize> dim = {(mwSize)(dim0), (mwSize)(1)};
    
    double *array = initAndAttachArray(matlabStruct, fieldName, dim);
    
    for(int i = 0; i < dim0; i++)
        array[i] = static_cast<double>(fieldData[i]);
}

template<typename T>
void writeMatlabField2(mxArray *matlabStruct, const char *fieldName,
                      const std::vector<T> fieldData, int dim0, int dim1,
                      std::vector<int> perm) {
    /**
     * @brief initialise matrix, attach to the field and write data
     * @param matlabStruct Pointer to the matlab structure
     * @param fieldName Name of the field to which the tensor will be attached
     * @param fieldData Data wich will be stored in the field
     * @param dim0 number of rows in the tensor
     * @param dim1 number of columns in the tensor
     */
    if(fieldData.size() != dim0*dim1)
        throw AmiException("Dimension mismatch when writing rdata->%s to matlab results",fieldName);
    
    if(perm.size() != 2)
        throw AmiException("Dimension mismatch when applying permutation!");
    
    std::vector<mwSize> dim = {(mwSize)(dim0), (mwSize)(dim1)};
    
    double *array = initAndAttachArray(matlabStruct, fieldName, reorder(dim,perm));
    
    std::vector<int> index = {0,0};
    /* transform rowmajor (c++) to colmajor (matlab) and apply permutation */
    for (index[0] = 0; index[0] < dim[0]; index[0]++) {
        for (index[1] = 0; index[1] < dim[1]; index[1]++) {
            array[index[perm[0]] + index[perm[1]]*dim[perm[0]]] =
                static_cast<double>(fieldData[index[0]*dim[1] + index[1]]);
        }
    }
}

template<typename T>
void writeMatlabField3(mxArray *matlabStruct, const char *fieldName,
                      const std::vector<T> fieldData, int dim0, int dim1,
                      int dim2, std::vector<int> perm) {
    /**
     * @brief initialise 3D tensor, attach to the field and write data
     * @param matlabStruct Pointer to the matlab structure
     * @param fieldName Name of the field to which the tensor will be attached
     * @param fieldData Data wich will be stored in the field
     * @param dim0 number of rows in the tensor
     * @param dim1 number of columns in the tensor
     * @param dim2 number of elements in the third dimension of the tensor
     */
    if(fieldData.size() != dim0*dim1*dim2)
        throw AmiException("Dimension mismatch when writing rdata->%s to matlab results",fieldName);
    
    if(perm.size() != 3)
        throw AmiException("Dimension mismatch when applying permutation!");
    
    std::vector<mwSize> dim = {(mwSize)(dim0), (mwSize)(dim1), (mwSize)(dim2)};
    
    double *array = initAndAttachArray(matlabStruct, fieldName, reorder(dim,perm));
    
    std::vector<int> index = {0,0,0};
    /* transform rowmajor (c++) to colmajor (matlab) and apply permutation */
    for (index[0] = 0; index[0] < dim[0]; index[0]++) {
        for (index[1] = 0; index[1] < dim[1]; index[1]++) {
            for (index[2] = 0; index[2] < dim[2]; index[2]++) {
                array[index[perm[0]] + (index[perm[1]] + index[perm[2]]*dim[perm[1]])*dim[perm[0]]] =
                    static_cast<double>(fieldData[(index[0]*dim[1] + index[1])*dim[2] + index[2]]);
            }
        }
    }
}

template<typename T>
void writeMatlabField4(mxArray *matlabStruct, const char *fieldName,
                      const std::vector<T> fieldData, int dim0, int dim1,
                      int dim2, int dim3, std::vector<int> perm) {
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
    if(fieldData.size() != dim0*dim1*dim2*dim3)
        throw AmiException("Dimension mismatch when writing rdata->%s to matlab results!",fieldName);
    
    if(perm.size() != 4)
        throw AmiException("Dimension mismatch when applying permutation!");
    
    std::vector<mwSize> dim = {(mwSize)(dim0), (mwSize)(dim1), (mwSize)(dim2), (mwSize)(dim3)};
    
    double *array = initAndAttachArray(matlabStruct, fieldName, reorder(dim,perm));
    
    std::vector<int> index = {0,0,0,0};
    /* transform rowmajor (c++) to colmajor (matlab) and apply permutation */
    for (index[0] = 0; index[0] < dim[0]; index[0]++) {
        for (index[1] = 0; index[1] < dim[1]; index[1]++) {
            for (index[2] = 0; index[2] < dim[2]; index[2]++) {
                for (index[3] = 0; index[3] < dim[3]; index[3]++) {
                    array[index[perm[0]] + (index[perm[1]] + (index[perm[2]] + index[perm[3]]*dim[perm[2]])*dim[perm[1]])*dim[perm[0]]] =
                        static_cast<double>(fieldData[((index[0]*dim[1] + index[1])*dim[2] + index[2])*dim[3] + index[3]]);
                }
            }
        }
    }
}

double *initAndAttachArray(mxArray *matlabStruct, const char *fieldName, std::vector<mwSize> dim) {
    if(!mxIsStruct(matlabStruct))
        throw AmiException("Passing non-struct mxArray to initAndAttachArray!",fieldName);
    
    int fieldNumber = mxGetFieldNumber(matlabStruct, fieldName);
    if(fieldNumber<0)
        throw AmiException("Trying to access non-existent field '%s'!",fieldName);

    mxArray *array = mxCreateNumericArray(dim.size(), dim.data(), mxDOUBLE_CLASS, mxREAL);
    mxSetFieldByNumber(matlabStruct, 0, fieldNumber, array);
    return(mxGetPr(array));
}

void checkFieldNames(const char **fieldNames,const int fieldCount) {
    for (int ifield = 0; ifield<fieldCount; ifield++) {
        if(!fieldNames[ifield])
            throw AmiException("Incorrect field name allocation, number of fields is smaller than fieldCount!");
    }
}

template<typename T>
std::vector<T> reorder(const std::vector<T> input, const std::vector<int> order) {
    if(order.size() != input.size())
        throw AmiException("Input dimension mismatch!");
    std::vector<T> reordered;
    reordered.resize(input.size());
    for(int i = 0; i < input.size(); i++)
        reordered[i] = input[order[i]];
    return(reordered);
}


} // namespace amici
