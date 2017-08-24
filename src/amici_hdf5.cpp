#include "include/amici_hdf5.h"
#include "include/amici_interface_cpp.h"
#include "include/amici_model.h"

#include <cassert>
#include <cstring>
#ifdef AMI_HDF5_H_DEBUG
#ifndef __APPLE__
#include <cexecinfo>
#else
#include <execinfo.h>
#endif
#endif
#include <unistd.h>

UserData *AMI_HDF5_readSimulationUserDataFromFileName(const char *fileName,
                                                      const char *datasetPath,
                                                      Model *model) {

    hid_t file_id = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);

    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileObject(
        file_id, datasetPath, model);

    H5Fclose(file_id);

    return (udata);
}

UserData *AMI_HDF5_readSimulationUserDataFromFileObject(hid_t fileId,
                                                        const char *datasetPath,
                                                        Model *model) {
    assert(fileId > 0);

    UserData *udata = new UserData();

    if (udata == NULL)
        return (NULL);

    udata->atol =
        AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "atol");
    udata->rtol =
        AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "rtol");
    udata->maxsteps =
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "maxsteps");
    udata->tstart =
        AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "tstart");
    udata->lmm = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "lmm");
    udata->iter = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "iter");
    udata->linsol =
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "linsol");
    udata->stldet =
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "stldet");
    udata->interpType =
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "interpType");
    udata->ism = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "ism");
    udata->sensi_meth = (AMICI_sensi_meth)AMI_HDF5_getIntScalarAttribute(
        fileId, datasetPath, "sensi_meth");
    udata->sensi = (AMICI_sensi_order)AMI_HDF5_getIntScalarAttribute(
        fileId, datasetPath, "sensi");
    udata->nmaxevent =
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "nmaxevent");
    udata->ordering =
        AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "ordering");

    hsize_t length;
    int status = 0;

    status += AMI_HDF5_getDoubleArrayAttribute(
        fileId, datasetPath, "qpositivex", &udata->qpositivex, &length);
    if (length != (unsigned)model->nx)
        return NULL;

    if (AMI_HDF5_attributeExists(fileId, datasetPath, "theta")) {
        status += AMI_HDF5_getDoubleArrayAttribute(fileId, datasetPath, "theta",
                                                   &udata->p, &length);
        if ((unsigned)model->np != length)
            return NULL;
    }

    if (AMI_HDF5_attributeExists(fileId, datasetPath, "kappa")) {
        status += AMI_HDF5_getDoubleArrayAttribute(fileId, datasetPath, "kappa",
                                                   &udata->k, &length);
        if (length != (unsigned)model->nk)
            return NULL;
    }

    if (AMI_HDF5_attributeExists(fileId, datasetPath, "ts")) {
        status += AMI_HDF5_getDoubleArrayAttribute(fileId, datasetPath, "ts",
                                                   &udata->ts, &length);
        udata->nt = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "nt");
        if (length != (unsigned)udata->nt || status > 0)
            return NULL;
    }

    if (AMI_HDF5_attributeExists(fileId, datasetPath, "pscale")) {
        udata->pscale = (AMICI_parameter_scaling)AMI_HDF5_getIntScalarAttribute(
            fileId, datasetPath, "pscale");
    }

    // parameter selection and reordering for sensitivities (matlab: fifth
    // argument)
    AMI_HDF5_getIntArrayAttribute(fileId, datasetPath, "sens_ind",
                                  &udata->plist, &length);
    assert(udata->nplist <= model->np);

    if (length > 0) {
        udata->nplist = length;
        // TODO: currently base 1 indices are written
        for (int i = 0; i < udata->nplist; ++i)
            udata->plist[i] -= 1;
    } else {
        udata->nplist = model->np;
        udata->plist = new int[udata->nplist];
        for (int i = 0; i < model->np; ++i)
            udata->plist[i] = i;
    }

    /* Options ; matlab: fourth argument   */
    // user-provided sensitivity initialisation. this should be a matrix of
    // dimension [#states x #parameters] default is sensitivity initialisation
    // based on the derivative of the state initialisation
    udata->x0data = NULL;
    udata->sx0data = NULL;

    /* pbarm parameterscales ; matlab: sixth argument*/
    udata->pbar = NULL;

    // xscale, matlab: seventh argument
    udata->xbar = NULL;

    return udata;
}

ExpData *AMI_HDF5_readSimulationExpData(const char *hdffile, UserData *udata,
                                        const char *dataObject, Model *model) {

    ExpData *edata = new ExpData(udata, model);

    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t m, n;

    if (H5Lexists(file_id, dataObject, 0)) {
        double *tmp_data;
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Y", &tmp_data,
                                           &m, &n);
        // if this is rank 1, n and m can be swapped
        if (n == 1) {
            assert(n == (unsigned)udata->nt || n == (unsigned)model->nytrue);
            assert(m == (unsigned)model->nytrue || m == (unsigned)udata->nt);
            assert(m * n == (unsigned)model->nytrue * udata->nt);
        } else {
            assert(n == (unsigned)udata->nt);
            assert(m == (unsigned)model->nytrue);
        }
        memcpy(edata->my, tmp_data, udata->nt * model->nytrue * sizeof(double));
        delete[] tmp_data;

        AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Sigma_Y",
                                           &tmp_data, &m, &n);
        if (n == 1) {
            assert(n == (unsigned)udata->nt || n == (unsigned)model->nytrue);
            assert(m == (unsigned)model->nytrue || m == (unsigned)udata->nt);
            assert(m * n == (unsigned)model->nytrue * udata->nt);
        } else {
            assert(n == (unsigned)udata->nt);
            assert(m == (unsigned)model->nytrue);
        }
        memcpy(edata->sigmay, tmp_data,
               udata->nt * model->nytrue * sizeof(double));
        delete[] tmp_data;

        if (model->nz) {
            AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Z",
                                               &tmp_data, &m, &n);
            if (n == 1) {
                assert(n == (unsigned)udata->nmaxevent ||
                       n == (unsigned)model->nztrue);
                assert(m == (unsigned)model->nztrue ||
                       m == (unsigned)udata->nmaxevent);
                assert(m * n == (unsigned)model->nytrue * udata->nmaxevent);
            } else {
                assert(n == (unsigned)udata->nmaxevent);
                assert(m == (unsigned)model->nztrue);
            }
            memcpy(edata->mz, tmp_data,
                   udata->nmaxevent * model->nztrue * sizeof(double));
            delete[] tmp_data;

            AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Sigma_Z",
                                               &tmp_data, &m, &n);
            if (n == 1) {
                assert(n == (unsigned)udata->nmaxevent ||
                       n == (unsigned)model->nztrue);
                assert(m == (unsigned)model->nztrue ||
                       m == (unsigned)udata->nmaxevent);
                assert(m * n == (unsigned)model->nytrue * udata->nmaxevent);
            } else {
                assert(n == (unsigned)udata->nmaxevent);
                assert(m == (unsigned)model->nztrue);
            }
            memcpy(edata->sigmaz, tmp_data,
                   udata->nmaxevent * model->nztrue * sizeof(double));
            delete[] tmp_data;
        }
    }
    H5Fclose(file_id);

    return edata;
}

void AMI_HDF5_writeReturnData(const ReturnData *rdata, const UserData *udata,
                              const char *hdffile, const char *datasetPath) {

    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDWR, H5P_DEFAULT);

    hid_t dataset;

    if (H5Lexists(file_id, datasetPath, H5P_DEFAULT)) {
        printf("INFO: result section already exists -- overwriting.\n");
        dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);
    } else {
        hsize_t dim[] = {1};
        hid_t dataspace = H5Screate_simple(1, dim, NULL);
        dataset = H5Dcreate(file_id, datasetPath, H5T_IEEE_F64LE, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (rdata->ts)
        H5LTset_attribute_double(file_id, datasetPath, "t", rdata->ts,
                                 rdata->nt);

    if (rdata->xdot)
        H5LTset_attribute_double(file_id, datasetPath, "xdot", rdata->xdot,
                                 rdata->nx);

    if (rdata->llh)
        H5LTset_attribute_double(file_id, datasetPath, "llh", rdata->llh, 1);

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
                                         const char *attributeName) {
    double doubleScalar;

    herr_t status = H5LTget_attribute_double(file_id, optionsObject,
                                             attributeName, &doubleScalar);
    assert(status >= 0);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %e\n", attributeName, doubleScalar);
#endif

    return doubleScalar;
}

int AMI_HDF5_getIntScalarAttribute(hid_t file_id, const char *optionsObject,
                                   const char *attributeName) {
    int intScalar;

    H5LTget_attribute_int(file_id, optionsObject, attributeName, &intScalar);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d\n", attributeName, intScalar);
#endif

    return intScalar;
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

void AMI_HDF5_getDoubleArrayAttribute2D(hid_t file_id,
                                        const char *optionsObject,
                                        const char *attributeName,
                                        double **destination, hsize_t *m,
                                        hsize_t *n) {
    int rank;
    H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
    assert(rank <= 2);

    hsize_t dims[2];
    H5T_class_t type_class;
    size_t type_size;
    H5LTget_attribute_info(file_id, optionsObject, attributeName, dims,
                           &type_class, &type_size);

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
}

void AMI_HDF5_getDoubleArrayAttribute3D(hid_t file_id,
                                        const char *optionsObject,
                                        const char *attributeName,
                                        double **destination, hsize_t *m,
                                        hsize_t *n, hsize_t *o) {
    int rank;
    H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
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
