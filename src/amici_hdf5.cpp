#include "include/amici_hdf5.h"
#include "include/amici.h"
#include "include/amici_model.h"
#include "include/edata.h"
#include "include/rdata.h"
#include <amici_solver.h>
#include <amici_exception.h>

#include <cassert>
#include <cstring>
#include <algorithm>
#ifdef AMI_HDF5_H_DEBUG
#ifndef __APPLE__
#include <cexecinfo>
#else
#include <execinfo.h>
#endif
#endif
#include <unistd.h>

void getModelDims(int *nx, int *nk, int *np);

namespace amici {
ExpData *AMI_HDF5_readSimulationExpData(const char *hdffile,
                                        const char *dataObject, Model *model) {

    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    ExpData *edata = NULL;

    hsize_t m, n;

    if (H5Lexists(file_id, dataObject, 0)) {

        edata = new ExpData(model);

        double *tmp_data;
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Y", &tmp_data,
                                           &m, &n);
        // if this is rank 1, n and m can be swapped
        if (n == 1) {
            assert(n == (unsigned)model->nt() || n == (unsigned)model->nytrue);
            assert(m == (unsigned)model->nytrue || m == (unsigned)model->nt());
            assert(m * n == (unsigned)model->nytrue * model->nt());
        } else {
            assert(n == (unsigned)model->nt());
            assert(m == (unsigned)model->nytrue);
        }
        edata->setObservedData(tmp_data);
        delete[] tmp_data;

        AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Sigma_Y",
                                           &tmp_data, &m, &n);
        if (n == 1) {
            assert(n == (unsigned)model->nt() || n == (unsigned)model->nytrue);
            assert(m == (unsigned)model->nytrue || m == (unsigned)model->nt());
            assert(m * n == (unsigned)model->nytrue * model->nt());
        } else {
            assert(n == (unsigned)model->nt());
            assert(m == (unsigned)model->nytrue);
        }
        edata->setObservedDataStdDev(tmp_data);
        delete[] tmp_data;

        if (model->nz) {
            AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Z",
                                               &tmp_data, &m, &n);
            if (n == 1) {
                assert(n == (unsigned)model->nMaxEvent() ||
                       n == (unsigned)model->nztrue);
                assert(m == (unsigned)model->nztrue ||
                       m == (unsigned)model->nMaxEvent());
                assert(m * n == (unsigned)model->nytrue * model->nMaxEvent());
            } else {
                assert(n == (unsigned)model->nMaxEvent());
                assert(m == (unsigned)model->nztrue);
            }
            edata->setObservedEvents(tmp_data);
            delete[] tmp_data;

            AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Sigma_Z",
                                               &tmp_data, &m, &n);
            if (n == 1) {
                assert(n == (unsigned)model->nMaxEvent() ||
                       n == (unsigned)model->nztrue);
                assert(m == (unsigned)model->nztrue ||
                       m == (unsigned)model->nMaxEvent());
                assert(m * n == (unsigned)model->nytrue * model->nMaxEvent());
            } else {
                assert(n == (unsigned)model->nMaxEvent());
                assert(m == (unsigned)model->nztrue);
            }
            edata->setObservedEventsStdDev(tmp_data);
            delete[] tmp_data;
        }
    }
    H5Fclose(file_id);

    return edata;
}

void AMI_HDF5_writeReturnData(const ReturnData *rdata,
                              const char *hdffile, const char *datasetPath) {

    hid_t file_id = H5Fcreate(hdffile, H5F_ACC_TRUNC , H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dim[] = {1};
    hid_t dataspace = H5Screate_simple(1, dim, NULL);

    auto plist_id  = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(plist_id, true);

    hid_t dataset = H5Dcreate(file_id, datasetPath, H5T_IEEE_F64LE, dataspace,
                              plist_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(plist_id);

    if (rdata->ts)
        H5LTset_attribute_double(file_id, datasetPath, "t", rdata->ts,
                                 rdata->nt);

    if (rdata->xdot)
        H5LTset_attribute_double(file_id, datasetPath, "xdot", rdata->xdot,
                                 rdata->nx);

    if (rdata->llh)
        H5LTset_attribute_double(file_id, datasetPath, "llh", rdata->llh, 1);
    
    if (rdata->status)
        H5LTset_attribute_double(file_id, datasetPath, "status", rdata->status, 1);

    if (rdata->sllh)
        H5LTset_attribute_double(file_id, datasetPath, "sllh", rdata->sllh,
                                 rdata->nplist);

    // are double, but should write as int:
    if (rdata->numsteps)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "numsteps",
                                           rdata->numsteps, rdata->nt);

    if (rdata->numrhsevals)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "numrhsevals",
                                           rdata->numrhsevals, rdata->nt);

    if (rdata->numerrtestfails)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath,
                                           "numerrtestfails",
                                           rdata->numerrtestfails, rdata->nt);

    if (rdata->numnonlinsolvconvfails)
        AMI_HDF5_setAttributeIntFromDouble(
            file_id, datasetPath, "numnonlinsolvconvfails",
            rdata->numnonlinsolvconvfails, rdata->nt);

    if (rdata->order)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "order",
                                           rdata->order, rdata->nt);

    if (rdata->numstepsB)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "numstepsB",
                                           rdata->numstepsB, rdata->nt);

    if (rdata->numrhsevalsB)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "numrhsevalsB",
                                           rdata->numrhsevalsB, rdata->nt);

    if (rdata->numerrtestfailsB)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath,
                                           "numerrtestfailsB",
                                           rdata->numerrtestfailsB, rdata->nt);

    if (rdata->numnonlinsolvconvfailsB)
        AMI_HDF5_setAttributeIntFromDouble(
            file_id, datasetPath, "numnonlinsolvconvfailsB",
            rdata->numnonlinsolvconvfailsB, rdata->nt);

    if (rdata->J)
        AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "J", rdata->J,
                                                 rdata->nx, rdata->nx);

    if (rdata->x)
        AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "x", rdata->x,
                                                 rdata->nt, rdata->nx);

    if (rdata->y)
        AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "y", rdata->y,
                                                 rdata->nt, rdata->ny);

    if (rdata->z)
        AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "z", rdata->z,
                                                 rdata->nmaxevent, rdata->nz);

    if (rdata->rz)
        AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "rz", rdata->rz,
                                                 rdata->nmaxevent, rdata->nz);

    if (rdata->sigmay)
        AMI_HDF5_createAndWriteDouble2DAttribute(
            dataset, "sigmay", rdata->sigmay, rdata->nt, rdata->ny);

    if (rdata->sigmaz)
        AMI_HDF5_createAndWriteDouble2DAttribute(
            dataset, "sigmaz", rdata->sigmaz, rdata->nt, rdata->nz);

    if (rdata->s2llh)
        AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "s2llh", rdata->s2llh,
                                                 rdata->nJ, rdata->nplist);

    if (rdata->sx)
        AMI_HDF5_createAndWriteDouble3DAttribute(
            dataset, "sx", rdata->sx, rdata->nt, rdata->nx, rdata->nplist);

    if (rdata->sy)
        AMI_HDF5_createAndWriteDouble3DAttribute(
            dataset, "sy", rdata->sy, rdata->nt, rdata->ny, rdata->nplist);

    if (rdata->ssigmay)
        AMI_HDF5_createAndWriteDouble3DAttribute(dataset, "ssigmay",
                                                 rdata->ssigmay, rdata->nt,
                                                 rdata->ny, rdata->nplist);

    if (rdata->sz)
        AMI_HDF5_createAndWriteDouble3DAttribute(dataset, "sz", rdata->sz,
                                                 rdata->nmaxevent, rdata->nz,
                                                 rdata->nplist);

    if (rdata->srz)
        AMI_HDF5_createAndWriteDouble3DAttribute(dataset, "srz", rdata->srz,
                                                 rdata->nmaxevent, rdata->nz,
                                                 rdata->nplist);

    if (rdata->ssigmaz)
        AMI_HDF5_createAndWriteDouble3DAttribute(
            dataset, "ssigmaz", rdata->ssigmaz, rdata->nmaxevent, rdata->nz,
            rdata->nplist);

    H5Fclose(file_id);
}

double AMI_HDF5_getDoubleScalarAttribute(hid_t file_id,
                                         const char *optionsObject,
                                         const char *attributeName, double *attributeValue) {
    herr_t status = H5LTget_attribute_double(file_id, optionsObject,
                                             attributeName, attributeValue);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %e\n", attributeName, *attributeValue);
#endif

    if(status < 0)
        warnMsgIdAndTxt("", "Attribute %s not found for object %s.", attributeName, optionsObject);

    return status >= 0;
}

int AMI_HDF5_getIntScalarAttribute(hid_t file_id, const char *optionsObject,
                                   const char *attributeName, int *attributeValue) {
    herr_t status = H5LTget_attribute_int(file_id, optionsObject, attributeName, attributeValue);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d\n", attributeName, *attributeValue);
#endif

    if(status < 0)
        warnMsgIdAndTxt("", "Attribute %s not found for object %s.", attributeName, optionsObject);

    return status >= 0;
}

int AMI_HDF5_getDoubleArrayAttribute(hid_t file_id, const char *optionsObject,
                                     const char *attributeName,
                                     double **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;
    herr_t status;
    *length = 0; // TODO: check why not set when reading scalar
    status = H5LTget_attribute_info(file_id, optionsObject, attributeName,
                                    length, &type_class, &type_size);

    if (status < 0) {
        fprintf(stderr, "Error in getDoubleArrayAttribute: Cannot read "
                        "attribute '%s' of '%s'\n",
                attributeName, optionsObject);

#ifdef AMI_HDF5_H_DEBUG
        void *array[10];
        size_t size;
        size = backtrace(array, 10);
        backtrace_symbols_fd(array, size, STDERR_FILENO);
#endif
        return status;
    }

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld: ", attributeName, *length);
#endif
    if(*length == 0) {
        *destination = NULL;
        return AMICI_SUCCESS;
    }

    *destination = new double[*length]; // vs. type_size
    status = H5LTget_attribute_double(file_id, optionsObject, attributeName,
                                      *destination);
    if (status < 0)
        fprintf(stderr, "Error in getDoubleArrayAttribute: Cannot read "
                        "attribute '%s' of '%s'\n",
                attributeName, optionsObject);

#ifdef AMI_HDF5_H_DEBUG
    // printfArray(*destination, *length, "%e ");
    printf("\n");
#endif
    return status < 0;
}

int AMI_HDF5_getDoubleArrayAttribute2D(hid_t file_id,
                                        const char *optionsObject,
                                        const char *attributeName,
                                        double **destination, hsize_t *m,
                                        hsize_t *n) {
    int rank;
    herr_t status = H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
    if(status < 0)
        return status;
    if(rank == 0) { // empty
        *destination = NULL;
        *m = *n = 0;
        return AMICI_SUCCESS;
    }

    assert(rank <= 2);
    hsize_t dims[2];
    H5T_class_t type_class;
    size_t type_size;
    status = H5LTget_attribute_info(file_id, optionsObject, attributeName, dims,
                           &type_class, &type_size);
    if(status < 0)
        return status;

    if (rank == 1) {
        *m = dims[0];
        *n = 1;
        AMI_HDF5_getDoubleArrayAttribute(file_id, optionsObject, attributeName,
                                         destination, m);
    } else {
#ifdef AMI_HDF5_H_DEBUG
        printf("%s: %lld x %lld: ", attributeName, dims[0], dims[1]);
#endif
        *m = dims[0];
        *n = dims[1];

        *destination = new double[(*m) * (*n)];
        H5LTget_attribute_double(file_id, optionsObject, attributeName,
                                 *destination);

#ifdef AMI_HDF5_H_DEBUG
        // printfArray(*destination, (*m) * (*n), "%e ");
        printf("\n");
#endif
    }
    return AMICI_SUCCESS;
}

int AMI_HDF5_getDoubleArrayAttribute3D(hid_t file_id, const char *optionsObject,
                                       const char *attributeName,
                                       double **destination, hsize_t *m,
                                       hsize_t *n, hsize_t *o) {
    if (!AMI_HDF5_attributeExists(file_id, optionsObject, attributeName))
        return 1;

    int rank;
    herr_t status =
        H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
    if (status < 0)
        return status;

    assert(rank == 3);

    hsize_t dims[3];
    H5T_class_t type_class;
    size_t type_size;
    H5LTget_attribute_info(file_id, optionsObject, attributeName, dims,
                           &type_class, &type_size);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld x %lld x %lld: ", attributeName, dims[0], dims[1],
           dims[2]);
#endif
    *m = dims[0];
    *n = dims[1];
    *o = dims[2];

    *destination = new double[(*m) * (*n) * (*o)];
    H5LTget_attribute_double(file_id, optionsObject, attributeName,
                             *destination);

#ifdef AMI_HDF5_H_DEBUG
    // printfArray(*destination, (*m) * (*n) * (*o), "%e ");
    printf("\n");
#endif

    return AMICI_SUCCESS;
}

void AMI_HDF5_getDoubleArrayAttribute4D(hid_t file_id,
                                        const char *optionsObject,
                                        const char *attributeName,
                                        double **destination, hsize_t *m,
                                        hsize_t *n, hsize_t *o, hsize_t *pp) {
    int rank;
    H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
    assert(rank == 4);

    hsize_t dims[4];
    H5T_class_t type_class;
    size_t type_size;
    H5LTget_attribute_info(file_id, optionsObject, attributeName, dims,
                           &type_class, &type_size);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld x %lld x %lld x %lld: ", attributeName, dims[0], dims[1],
           dims[2], dims[3]);
#endif
    *m = dims[0];
    *n = dims[1];
    *o = dims[2];
    *pp = dims[3];

    *destination = new double[(*m) * (*n) * (*o) * (*pp)];
    H5LTget_attribute_double(file_id, optionsObject, attributeName,
                             *destination);

#ifdef AMI_HDF5_H_DEBUG
    // printfArray(*destination, (*m) * (*n) * (*o) * (*pp), "%e ");
    printf("\n");
#endif
}

void AMI_HDF5_getIntArrayAttribute(hid_t file_id, const char *optionsObject,
                                   const char *attributeName, int **destination,
                                   hsize_t *length) {
    *length = 0;

    H5T_class_t type_class;
    size_t type_size;

    H5LTget_attribute_info(file_id, optionsObject, attributeName, length,
                           &type_class, &type_size);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld: ", attributeName, *length);
#endif
    if (*length > 0) {
        *destination = new int[*length];
        H5LTget_attribute_int(file_id, optionsObject, attributeName,
                              *destination);
    } else {
        *destination = NULL;
    }
}

// TODO: option for overwrite
herr_t AMI_HDF5_createAndWriteDouble2DAttribute(hid_t dataset,
                                                const char *attributeName,
                                                const double *buffer, hsize_t m,
                                                hsize_t n) {
    const hsize_t adims[] = {m, n};

    if (H5Aexists(dataset, attributeName) > 0) {
        H5Adelete(dataset, attributeName);
    }

    hid_t space = H5Screate_simple(2, adims, NULL);
    hid_t attr = H5Acreate2(dataset, attributeName, H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr, H5T_NATIVE_DOUBLE, buffer);

    return status;
}

herr_t AMI_HDF5_createAndWriteDouble3DAttribute(hid_t dataset,
                                                const char *attributeName,
                                                const double *buffer, hsize_t m,
                                                hsize_t n, hsize_t o) {
    const hsize_t adims[] = {m, n, o};

    if (H5Aexists(dataset, attributeName) > 0) {
        H5Adelete(dataset, attributeName);
    }

    hid_t space = H5Screate_simple(3, adims, NULL);
    hid_t attr = H5Acreate2(dataset, attributeName, H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr, H5T_NATIVE_DOUBLE, buffer);

    return status;
}

void AMI_HDF5_setAttributeIntFromDouble(hid_t file_id, const char *obj_name,
                                        const char *attr_name,
                                        const double *bufferDouble,
                                        size_t size) {
    int intBuffer[size];
    for (int i = 0; (unsigned)i < size; ++i)
        intBuffer[i] = bufferDouble[i];

    H5LTset_attribute_int(file_id, obj_name, attr_name, intBuffer, size);
}

int AMI_HDF5_attributeExists(hid_t fileId, const char *datasetPath,
                             const char *attributeName) {
    if (H5Lexists(fileId, datasetPath, H5P_DEFAULT)) {
        hid_t dataset = H5Dopen2(fileId, datasetPath, H5P_DEFAULT);
        if (H5LTfind_attribute(dataset, attributeName))
            return 1;
    }

    return 0;
}

void readSolverSettingsFromHDF5(hid_t fileId, Solver &solver, const std::string &datasetPath_) {
    const char *datasetPath = datasetPath_.c_str();
    double dblBuffer = 0;
    int intBuffer = 0;

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "atol")) {
        AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "atol", &dblBuffer);
        solver.setAbsoluteTolerance(dblBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "rtol")) {
        AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "rtol", &dblBuffer);
        solver.setRelativeTolerance(dblBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "maxsteps")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "maxsteps", &intBuffer);
        solver.setMaxSteps(intBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "lmm")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "lmm", &intBuffer);
        solver.setLinearMultistepMethod(static_cast<LinearMultistepMethod>(intBuffer));
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "iter")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "iter", &intBuffer);
        solver.setNonlinearSolverIteration(static_cast<NonlinearSolverIteration>(intBuffer));
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "stldet")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "stldet", &intBuffer);
        solver.setStabilityLimitFlag(intBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "ordering")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "ordering", &intBuffer);
        solver.setStateOrdering(static_cast<StateOrdering>(intBuffer));
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "interpType")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "interpType", &intBuffer);
        solver.setInterpolationType(static_cast<InterpolationType>(intBuffer));
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "sensi_meth")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "sensi_meth", &intBuffer);
        solver.setSensitivityMethod(static_cast<AMICI_sensi_meth>(intBuffer));
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "sensi")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "sensi", &intBuffer);
        solver.setSensitivityOrder(static_cast<AMICI_sensi_order>(intBuffer));
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "newton_maxsteps")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "newton_maxsteps", &intBuffer);
        solver.setNewtonMaxSteps(intBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "newton_preeq")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "newton_preeq", &intBuffer);
        solver.setNewtonPreequilibration(intBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "newton_precon")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "newton_precon", &intBuffer);
        solver.setNewtonPreconditioner(intBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "newton_maxlinsteps")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "newton_maxlinsteps", &intBuffer);
        solver.setNewtonMaxLinearSteps(intBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "linsol")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "linsol", &intBuffer);
        solver.setLinearSolver(static_cast<LinearSolver>(intBuffer));
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "ism")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "ism", &intBuffer);
        solver.setInternalSensitivityMethod(static_cast<InternalSensitivityMethod>(intBuffer));
    }
}

void readSolverSettingsFromHDF5(const std::string &hdffile, Solver &solver, const std::string &datasetPath) {
    hid_t file_id = H5Fopen(hdffile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    readSolverSettingsFromHDF5(file_id, solver, datasetPath);

    H5Fclose(file_id);
}

void readModelDataFromHDF5(const std::string &hdffile, Model &model, const std::string &datasetPath) {
    hid_t file_id = H5Fopen(hdffile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    readModelDataFromHDF5(file_id, model, datasetPath);

    H5Fclose(file_id);
}

void readModelDataFromHDF5(hid_t fileId, Model &model, const std::string &datasetPath_) {
    const char *datasetPath = datasetPath_.c_str();
    int intBuffer = 0;
    double dblBuffer = 0;

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "tstart")) {
        AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "tstart", &dblBuffer);
        model.setT0(dblBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "pscale")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "pscale", &intBuffer);
        model.setParameterScale(static_cast<AMICI_parameter_scaling>(intBuffer));
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "nmaxevent")) {
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "nmaxevent", &intBuffer);
        model.setNMaxEvent(intBuffer);
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "qpositivex")) {
        double *buffer = nullptr;
        hsize_t length0 = 0;
        // TODO double vs int?!
        AMI_HDF5_getDoubleArrayAttribute(
                    fileId, datasetPath, "qpositivex", &buffer, &length0);
        if (length0 == (unsigned) model.nx)
            model.setPositivityFlag(std::vector<int>(buffer, buffer + length0));
        else if(length0 != 0) { // currently not written correctly from matlab
            delete[] buffer;
            throw(AmiException("Failed reading qpositivex (%d).", length0));
        }
        delete[] buffer;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "theta")) {
        double *buffer = nullptr;
        hsize_t length0 = 0;
        AMI_HDF5_getDoubleArrayAttribute(
                    fileId, datasetPath, "theta", &buffer, &length0);
        if (length0 != (unsigned) model.np()) {
            delete[] buffer;
            throw(AmiException("Failed reading theta (%d).", length0));
        }
        model.setParameters(std::vector<double>(buffer, buffer + length0));
        delete[] buffer;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "kappa")) {
        double *buffer = nullptr;
        hsize_t length0 = 0;
        AMI_HDF5_getDoubleArrayAttribute(
                    fileId, datasetPath, "kappa", &buffer, &length0);
        if (length0 != (unsigned) model.nk()) {
            delete[] buffer;
            throw(AmiException("Failed reading kappa."));
        }
        model.setFixedParameters(std::vector<double>(buffer, buffer + length0));
        delete[] buffer;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "ts")) {
        double *buffer = nullptr;
        hsize_t length0 = 0;
        AMI_HDF5_getDoubleArrayAttribute(
                    fileId, datasetPath, "ts", &buffer, &length0);
        model.setTimepoints(std::vector<double>(buffer, buffer + length0));
        delete[] buffer;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "sens_ind")) {
        int *buffer = nullptr;
        hsize_t length0 = 0;
        AMI_HDF5_getIntArrayAttribute(
                    fileId, datasetPath, "sens_ind", &buffer, &length0);
        if (length0 > 0) {
            // currently base 1 indices are written
            for (int i = 0; (unsigned)i < length0; ++i) {
                buffer[i] -= 1;
            }
            model.setParameterList(std::vector<int>(buffer, buffer + length0));
        } else {
            model.requireSensitivitiesForAllParameters();
        }
        delete[] buffer;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "x0")) {
        double *buffer = nullptr;
        hsize_t length0 = 0;
        AMI_HDF5_getDoubleArrayAttribute(
                    fileId, datasetPath, "x0", &buffer, &length0);
        if (length0 == (unsigned) model.nx)
            model.setInitialStates(std::vector<double>(buffer, buffer + length0));
        delete[] buffer;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "sx0")) {
        double *buffer = nullptr;
        hsize_t length0 = 0;
        hsize_t length1 = 0;
        AMI_HDF5_getDoubleArrayAttribute2D(
                    fileId, datasetPath, "sx0", &buffer, &length0, &length1);
        if (length0 == (unsigned) model.nx && length1 == (unsigned) model.nplist())
            model.setInitialStateSensitivities(std::vector<double>(buffer, buffer + length0));
        delete[] buffer;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "pbar")) {
        double *buffer = nullptr;
        hsize_t length0 = 0;
        AMI_HDF5_getDoubleArrayAttribute(
                    fileId, datasetPath, "pbar", &buffer, &length0);
        if (length0 == (unsigned) model.nplist())
            model.setParameterScaling(std::vector<double>(buffer, buffer + length0));
        delete[] buffer;
    }
}

} // namespace amici
