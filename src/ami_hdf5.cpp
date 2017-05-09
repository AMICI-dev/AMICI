#include "include/ami_hdf5.h"
#include "include/amici_interface_cpp.h"

#include <assert.h>
#ifdef AMI_HDF5_H_DEBUG
#include <execinfo.h>
#endif
#include <unistd.h>

#include "wrapfunctions.h"
#include "include/symbolic_functions.h"

#include "include/rdata_accessors.h"
#include "include/edata_accessors.h"
#include "include/udata_accessors.h"


UserData *AMI_HDF5_readSimulationUserDataFromFileName(const char* fileName, const char* datasetPath) {

    hid_t file_id = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);

    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileObject(file_id, datasetPath);

    H5Fclose(file_id);

    return(udata);
}

UserData *AMI_HDF5_readSimulationUserDataFromFileObject(hid_t fileId, const char *datasetPath)
{
    UserData *udata = getDefaultUserData();
    if (udata == NULL)
        return(NULL);

    // Set hard-coded model dimensions
    init_modeldims(udata);

    atol       = AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "atol");
    rtol       = AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "rtol");
    maxsteps   = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "maxsteps");
    tstart     = AMI_HDF5_getDoubleScalarAttribute(fileId, datasetPath, "tstart");
    lmm        = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "lmm");
    iter       = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "iter");
    linsol     = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "linsol");
    stldet     = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "stldet");
    interpType = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "interpType");
    ism        = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "ism");
    sensi_meth = (AMI_sensi_meth) AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "sensi_meth");
    sensi      = (AMI_sensi_order) AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "sensi");
    nmaxevent  = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "nmaxevent");
    ordering   = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "ordering");

    hsize_t length;
    int status = 0;

    status += AMI_HDF5_getDoubleArrayAttribute(fileId, datasetPath, "qpositivex", &qpositivex, &length);

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "theta")) {
        status += AMI_HDF5_getDoubleArrayAttribute(fileId, datasetPath, "theta", &p, &length);
        if(np != length)
            return NULL;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "kappa")) {
        status += AMI_HDF5_getDoubleArrayAttribute(fileId, datasetPath, "kappa", &k, &length);
        if(length != nk)
            return NULL;
    }

    if(AMI_HDF5_attributeExists(fileId, datasetPath, "ts")) {
        status += AMI_HDF5_getDoubleArrayAttribute(fileId, datasetPath, "ts", &ts, &length);
        nt = AMI_HDF5_getIntScalarAttribute(fileId, datasetPath, "nt");
        if(length != nt || status > 0)
            return NULL;
    }

    // parameter selection and reordering for sensitivities (matlab: fifth argument)
    // For now, use all parameters
    nplist = np;
    plist = new int[np];
    for (int i = 0; i < np; i++)
        plist[i] = i;

    /* Options ; matlab: fourth argument   */
    z2event = new realtype[ne];
    for(int i = 0; i < ne; ++i)
        z2event[i] = i;

    idlist = new realtype[nplist]();
    for(int i = 0; i < nplist; ++i)
        idlist[i] = 0;

    //user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters] default is sensitivity initialisation based on the derivative of the state initialisation
    x0data = NULL;
    sx0data = NULL;

    /* pbarm parameterscales ; matlab: sixth argument*/
    pbar = NULL;

    // xscale, matlab: seventh argument
    xbar = NULL;

    h = 0;

    udata->am_dxdotdp = NULL;
    udata->am_dwdx = NULL;
    udata->am_dwdp = NULL;
    udata->am_M = NULL;
    udata->am_dfdx = NULL;
    udata->am_stau = NULL;

    processUserData(udata);

    return udata;
}


ExpData *AMI_HDF5_readSimulationExpData(const char* hdffile, UserData *udata, const char* dataObject) {
    ExpData *edata = new ExpData();
    if (edata == NULL) {
        return(NULL);
    }

    mz = NULL;
    zsigma = NULL;

    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t m, n;

    if(H5Lexists(file_id, dataObject, 0)) {
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Y", &my, &m, &n);
        assert(m * n == nt * ny); // TODO m, n separately

        AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Sigma_Y", &ysigma, &m, &n);
        assert(m * n == nt * ny);

        if(nz) {
            AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Z", &mz, &m, &n);
            assert(m * n == nt * nz);

            AMI_HDF5_getDoubleArrayAttribute2D(file_id, dataObject, "Sigma_Z", &zsigma, &m, &n);
            assert(m * n == nt * nz);
        }
    }
    H5Fclose(file_id);

    return(edata);
}

void AMI_HDF5_writeReturnData(const ReturnData *rdata, const UserData *udata, const char* hdffile, const char* datasetPath) {

    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDWR, H5P_DEFAULT);

    hid_t dataset;

    if(H5Lexists(file_id, datasetPath, H5P_DEFAULT)) {
        printf("INFO: result section already exists -- overwriting.\n");
        dataset = H5Dopen2(file_id, datasetPath, H5P_DEFAULT);
    } else {
        hsize_t dim [] = {1};
        hid_t dataspace = H5Screate_simple (1, dim, NULL);
        dataset = H5Dcreate(file_id, datasetPath, H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    H5LTset_attribute_double(file_id, datasetPath, "t", tsdata, nt);
    H5LTset_attribute_double(file_id, datasetPath, "xdot", xdotdata, nx);
    H5LTset_attribute_double(file_id, datasetPath, "llh", llhdata, 1);
    H5LTset_attribute_double(file_id, datasetPath, "sllh", sllhdata, np);

    // are double, but should write as int:
    AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "numsteps", numstepsdata, nt);
    AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "numrhsevals", numrhsevalsdata, nt);
    AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "order", orderdata, nt);

    if(numstepsSdata)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "numstepsS", numstepsSdata, nt);

    if(numrhsevalsSdata)
        AMI_HDF5_setAttributeIntFromDouble(file_id, datasetPath, "numrhsevalsS", numrhsevalsSdata, nt);

    AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "J", Jdata, nx, nx);
    AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "x", xdata, nt, nx);
    AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "y", ydata, nt, ny);
    AMI_HDF5_createAndWriteDouble2DAttribute(dataset, "sigmay", sigmaydata, nt, ny);

    if(sxdata)
        AMI_HDF5_createAndWriteDouble3DAttribute(dataset, "sx", sxdata, nt, nx, np);

    if(sydata)
        AMI_HDF5_createAndWriteDouble3DAttribute(dataset, "sy", sydata, nt, ny, np);

    if(ssigmaydata)
        AMI_HDF5_createAndWriteDouble3DAttribute(dataset, "ssigmay", ssigmaydata, nt, ny, np);

    H5Fclose(file_id);

    // NOT YET INCLUDED IN HDF5 output:
    //    /** parameter derivative of time derivative */
    //    double *am_dxdotdpdata;
    //    /** state derivative of observables */
    //    double *am_dydxdata;
    //    /** parameter derivative of observables */
    //    double *am_dydpdata;
    //    /** event output */
    //    double *am_zdata;
    //    /** event output sigma standard deviation */
    //    double *am_sigmazdata;
    //    /** parameter derivative of event output */
    //    double *am_szdata;
    //    /** parameter derivative of event output standard deviation */
    //    double *am_ssigmazdata;
    //    /** event trigger output */
    //    double *am_rzdata;
    //    /** parameter derivative of event trigger output */
    //    double *am_srzdata;
    //    /** second order parameter derivative of event trigger output */
    //    double *am_s2rzdata;
    //    /** chi2 value */
    //    double *am_chi2data;
    //    /** second order parameter derivative of likelihood */
    //    double *am_s2llhdata;
}


double AMI_HDF5_getDoubleScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName) {
    double doubleScalar;

    H5LTget_attribute_double(file_id, optionsObject, attributeName, &doubleScalar);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %e\n", attributeName, doubleScalar);
#endif

    return doubleScalar;
}

int AMI_HDF5_getIntScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName) {
    int intScalar;

    H5LTget_attribute_int(file_id, optionsObject, attributeName, &intScalar);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d\n", attributeName, intScalar);
#endif

    return intScalar;
}

int AMI_HDF5_getDoubleArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;
    herr_t status;
    *length = 0; // TODO: check why not set when reading scalar
    status = H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);

    if(status < 0) {
        fprintf(stderr, "Error in getDoubleArrayAttribute: Cannot read attribute '%s' of '%s'\n", attributeName, optionsObject);

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
    status = H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);
    if(status < 0)
        fprintf(stderr, "Error in getDoubleArrayAttribute: Cannot read attribute '%s' of '%s'\n", attributeName, optionsObject);

#ifdef AMI_HDF5_H_DEBUG
    printfArray(*destination, *length, "%e ");
    printf("\n");
#endif
    return status < 0;
}

void AMI_HDF5_getDoubleArrayAttribute2D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n) {
    int rank;
    H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
    assert(rank <= 2);

    hsize_t dims[2];
    H5T_class_t type_class;
    size_t type_size;
    H5LTget_attribute_info(file_id, optionsObject, attributeName, dims, &type_class, &type_size);

    if(rank == 1) {
        *m = 1;
        *n = dims[0];
        AMI_HDF5_getDoubleArrayAttribute(file_id, optionsObject, attributeName, destination, m);
    } else {
#ifdef AMI_HDF5_H_DEBUG
        printf("%s: %lld x %lld: ", attributeName, dims[0], dims[1]);
#endif
        *m = dims[0];
        *n = dims[1];

        *destination = new double[(*m) * (*n)];
        H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);

#ifdef AMI_HDF5_H_DEBUG
        printfArray(*destination, (*m) * (*n), "%e ");
        printf("\n");
#endif
    }
}

void AMI_HDF5_getDoubleArrayAttribute3D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n, hsize_t *o) {
    int rank;
    H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
    assert(rank == 3);

    hsize_t dims[3];
    H5T_class_t type_class;
    size_t type_size;
    H5LTget_attribute_info(file_id, optionsObject, attributeName, dims, &type_class, &type_size);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld x %lld x %lld: ", attributeName, dims[0], dims[1], dims[2]);
#endif
    *m = dims[0];
    *n = dims[1];
    *o = dims[2];

    *destination = new double[(*m) * (*n) * (*o)];
    H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);

#ifdef AMI_HDF5_H_DEBUG
    printfArray(*destination, (*m) * (*n) * (*o), "%e ");
    printf("\n");
#endif
}

void AMI_HDF5_getDoubleArrayAttribute4D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n, hsize_t *o, hsize_t *pp) {
    int rank;
    H5LTget_attribute_ndims(file_id, optionsObject, attributeName, &rank);
    assert(rank == 4);

    hsize_t dims[4];
    H5T_class_t type_class;
    size_t type_size;
    H5LTget_attribute_info(file_id, optionsObject, attributeName, dims, &type_class, &type_size);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld x %lld x %lld x %lld: ", attributeName, dims[0], dims[1], dims[2], dims[3]);
#endif
    *m = dims[0];
    *n = dims[1];
    *o = dims[2];
    *pp = dims[3];

    *destination = new double[(*m) * (*n) * (*o) * (*pp)];
    H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);

#ifdef AMI_HDF5_H_DEBUG
    printfArray(*destination, (*m) * (*n) * (*o) * (*pp), "%e ");
    printf("\n");
#endif
}



void AMI_HDF5_getIntArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;

    H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);

#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %lld: ", attributeName, *length);
#endif

    *destination = new int[*length];
    H5LTget_attribute_int(file_id, optionsObject, attributeName, *destination);
}

// TODO: option for overwrite
herr_t AMI_HDF5_createAndWriteDouble2DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n) {
    const hsize_t adims[] = {m, n};

    if(H5Aexists(dataset, attributeName) > 0) {
        H5Adelete(dataset, attributeName);
    }

    hid_t space = H5Screate_simple(2, adims, NULL);
    hid_t attr = H5Acreate2(dataset, attributeName, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr, H5T_NATIVE_DOUBLE, buffer);

    return status;
}

herr_t AMI_HDF5_createAndWriteDouble3DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n, hsize_t o) {
    const hsize_t adims[] = {m, n, o};

    if(H5Aexists(dataset, attributeName) > 0) {
        H5Adelete(dataset, attributeName);
    }

    hid_t space = H5Screate_simple(3, adims, NULL);
    hid_t attr = H5Acreate2(dataset, attributeName, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr, H5T_NATIVE_DOUBLE, buffer);

    return status;
}


void AMI_HDF5_setAttributeIntFromDouble(hid_t file_id, const char *obj_name, const char *attr_name, const double *bufferDouble, size_t size)
{
    int intBuffer[size];
    for(int i = 0; i < size; ++i)
        intBuffer[i] = bufferDouble[i];

    H5LTset_attribute_int(file_id, obj_name, attr_name, intBuffer, size);
}

int AMI_HDF5_attributeExists(hid_t fileId, const char *datasetPath, const char *attributeName) {
    if(H5Lexists(fileId, datasetPath, H5P_DEFAULT)) {
        hid_t dataset = H5Dopen2(fileId, datasetPath, H5P_DEFAULT);
        if(H5LTfind_attribute(dataset, attributeName))
            return 1;
    }

    return 0;
}
