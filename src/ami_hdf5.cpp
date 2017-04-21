#include <assert.h>
#ifdef AMI_HDF5_H_DEBUG
#include <execinfo.h>
#endif AMI_HDF5_H_DEBUG
#include <unistd.h>

#include "ami_hdf5.h"
#include "wrapfunctions.h"
#include "include/symbolic_functions.h"

#include "include/rdata_accessors.h"
#include "include/edata_accessors.h"
#include "include/udata_accessors.h"

void storeSimulation(const char* fileName, ReturnData *rdata) {
    hsize_t dims[] = {1};

    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    int file_id;
    hid_t dataset_id = H5Dcreate(file_id, "/simulations", H5T_STD_I32BE,
                            dataspace_id,
                            H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
    hid_t status;
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
}

UserData *readSimulationUserData(const char* fileName) {
    UserData *udata = new UserData();
    if (udata == NULL)
        return(NULL);

    init_modeldims(udata);

    hid_t file_id = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    //status = H5Aread (attr_id, mem_type_id, buf);

    const char* optionsObject = "/options";
    atol = getDoubleScalarAttribute(file_id, optionsObject, "atol");
    rtol = getDoubleScalarAttribute(file_id, optionsObject, "rtol");
    maxsteps = getIntScalarAttribute(file_id, optionsObject, "maxsteps");
    tstart = getDoubleScalarAttribute(file_id, optionsObject, "tstart");
    lmm = getIntScalarAttribute(file_id, optionsObject, "lmm");
    iter = getIntScalarAttribute(file_id, optionsObject, "iter");
    linsol = getIntScalarAttribute(file_id, optionsObject, "linsol");
    stldet = getIntScalarAttribute(file_id, optionsObject, "stldet");
    interpType = getIntScalarAttribute(file_id, optionsObject, "interpType");
    //lmmB = getDoubleScalarAttribute(file_id, optionsObject, "lmmB");
    //iterB = getDoubleScalarAttribute(file_id, optionsObject, "iterB");
    ism = getIntScalarAttribute(file_id, optionsObject, "ism");
    sensi_meth = getIntScalarAttribute(file_id, optionsObject, "sensi_meth");
    sensi = getIntScalarAttribute(file_id, optionsObject, "sensi");
    nmaxevent = getIntScalarAttribute(file_id, optionsObject, "nmaxevent");
    ordering = getIntScalarAttribute(file_id, optionsObject, "ordering");
    //ss = getDoubleScalarAttribute(file_id, optionsObject, "ss");

    hsize_t length;
    getDoubleArrayAttribute(file_id, optionsObject, "qpositivex", &qpositivex, &length);
    // getIntArrayAttribute(file_id, optionsObject, "sens_ind", sens_ind, &tmp);
    //    'sx0':  []

    const char* dataObject = "/data";
    getDoubleArrayAttribute(file_id, dataObject, "theta", &p, &length);
    np = length;

    getDoubleArrayAttribute(file_id, dataObject, "kappa", &k, &length);
    assert(length == nk);
    getDoubleArrayAttribute(file_id, dataObject, "ts", &ts, &length);
    nt = getIntScalarAttribute(file_id, dataObject, "nt");
    assert(length == nt);
    assert(ny == getIntScalarAttribute(file_id, dataObject, "ny"));
    assert(nz == getIntScalarAttribute(file_id, dataObject, "nz"));


    /* plist, matlab: fifth argument */
    // parameter ordering
    plist = new int[np]();
    for (int i = 0; i < np; i++) {
        plist[i] = i;
    }

    /* Options ; matlab: fourth argument   */
    z2event = new realtype[ne]();
    for(int i = 0; i < ne; ++i)
        z2event[i] = i;

    idlist = new realtype[np]();
    for(int i = 0; i < np; ++i)
        idlist[i] = 0;
    //user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters] default is sensitivity initialisation based on the derivative of the state initialisation
    b_sx0 = FALSE;
    sx0data = 0;

    /* pbarm parameterscales ; matlab: sixth argument*/
    pbar = new realtype[np]();
    ones(pbar, np);

    //    /* xscale, matlab: seventh argument */
    //    xbar = mxGetPr(prhs[6]);
    xbar = 0; // xscale

//    double *am_qpositivex;
//    /** parameter reordering */
//    int    *am_plist;
//    /** number of parameters */
//    int    am_np;
//    /** number of observables */
//    int    am_ny;
//    /** number of observables in the unaugmented system */
//    int    am_nytrue;
//    /** number of states */
//    int    am_nx;
//    /** number of event outputs */
//    int    am_nz;
//    /** number of event outputs in the unaugmented system */
//    int    am_nztrue;
//    /** number of events */
//    int    am_ne;
//    /** number of timepoints */
//    int    am_nt;
//    /** number of common expressions */
//    int    am_nw;
//    /** number of derivatives of common expressions wrt x */
//    int    am_ndwdx;
//    /** number of derivatives of common expressions wrt p */
//    int    am_ndwdp;
//    /** number of nonzero entries in jacobian */
//    int    am_nnz;

//    /** scaling of parameters */
//    double *am_pbar;
//    /** scaling of states */
//    double *am_xbar;
//    /** flag array for DAE equations */
//    double *am_idlist;

//    /** upper bandwith of the jacobian */
//    int am_ubw;
//    /** lower bandwith of the jacobian */
//    int am_lbw;

//    /** flag for sensitivity initialisation */
//    /*!
//     * flag which determines whether analytic sensitivities initialisation or provided initialisation should be used
//     */
//    booleantype am_bsx0;

//    /** sensitivity initialisation */
//    double *am_sx0data;

//    /** index indicating to which event an event output belongs */
//    double *am_z2event;

//    /** flag indicating whether a certain heaviside function should be active or not */
//    double *am_h;

//    /** tempory storage of Jacobian data across functions */
//    SlsMat am_J;
//    /** tempory storage of dxdotdp data across functions */
//    realtype *am_dxdotdp;
//    /** tempory storage of w data across functions */
//    realtype *am_w;
//    /** tempory storage of dwdx data across functions */
//    realtype *am_dwdx;
//    /** tempory storage of dwdp data across functions */
//    realtype *am_dwdp;
//    /** tempory storage of M data across functions */
//    realtype *am_M;
//    /** tempory storage of dfdx data across functions */
//    realtype *am_dfdx;
//    /** tempory storage of stau data across functions */
//    realtype *am_stau;

//    /** flag indicating whether a NaN in dxdotdp has been reported */
//    booleantype am_nan_dxdotdp;
//    /** flag indicating whether a NaN in J has been reported */
//    booleantype am_nan_J;
//    /** flag indicating whether a NaN in JSparse has been reported */
//    booleantype am_nan_JSparse;
//    /** flag indicating whether a NaN in xdot has been reported */
//    booleantype am_nan_xdot;
//    /** flag indicating whether a NaN in xBdot has been reported */
//    booleantype am_nan_xBdot;
//    /** flag indicating whether a NaN in qBdot has been reported */
//    booleantype am_nan_qBdot;

    H5Fclose(file_id);
    processUserData(udata);

    return(udata);
}

ExpData *readSimulationExpData(const char* hdffile, UserData *udata) {
    ExpData *edata = new ExpData();
    if (edata == NULL) {
        return(NULL);
    }

    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t m, n;
    const char* dataObject = "/data";
    getDoubleArrayAttribute2D(file_id, dataObject, "Y", &my, &m, &n);
    assert(m * n == nt * ny); // TODO m, n separately
    getDoubleArrayAttribute2D(file_id, dataObject, "Sigma_Y", &ysigma, &m, &n);
    assert(m * n == nt * ny);
    if(nz) {
        getDoubleArrayAttribute2D(file_id, dataObject, "Z", &mz, &m, &n);
        assert(m * n == nt * nz);
        getDoubleArrayAttribute2D(file_id, dataObject, "Sigma_Z", &zsigma, &m, &n);
        assert(m * n == nt * nz);
    } else {
        mz = 0;
        zsigma = 0;
    }
    H5Fclose(file_id);

    b_expdata = TRUE;

    return(edata);
}

void writeReturnData(const char* hdffile, ReturnData *rdata, UserData *udata) {
    hid_t file_id = H5Fopen(hdffile, H5F_ACC_RDWR, H5P_DEFAULT);

    const char* solutionsObject = "/solutions";

    hid_t dataset;

    if(H5Lexists(file_id, solutionsObject, H5P_DEFAULT)) {
        printf("INFO: 'solutions' section already exists. overwriting.\n");
         dataset = H5Dopen2(file_id, solutionsObject, H5P_DEFAULT);
    } else {
        hsize_t dim [] = {1};
        hid_t dataspace = H5Screate_simple (1, dim, NULL);
        dataset = H5Dcreate(file_id, solutionsObject, H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    H5LTset_attribute_double(file_id, solutionsObject, "t", tsdata, nt);
    H5LTset_attribute_double(file_id, solutionsObject, "xdot", xdotdata, nx);
    H5LTset_attribute_double(file_id, solutionsObject, "llh", llhdata, 1);
    H5LTset_attribute_double(file_id, solutionsObject, "sllh", sllhdata, np);

    // are double, but should write as int:
    setAttributeIntFromDouble(file_id, solutionsObject, "numsteps", numstepsdata, nt);
    setAttributeIntFromDouble(file_id, solutionsObject, "numrhsevals", numrhsevalsdata, nt);
    setAttributeIntFromDouble(file_id, solutionsObject, "order", orderdata, nt);
    if(numstepsSdata)
        setAttributeIntFromDouble(file_id, solutionsObject, "numstepsS", numstepsSdata, nt);
    if(numrhsevalsSdata)
        setAttributeIntFromDouble(file_id, solutionsObject, "numrhsevalsS", numrhsevalsSdata, nt);

    createAndWriteDouble2DAttribute(dataset, "J", Jdata, nx, nx);
    createAndWriteDouble2DAttribute(dataset, "x", xdata, nt, nx);
    createAndWriteDouble2DAttribute(dataset, "y", ydata, nt, ny);
    createAndWriteDouble2DAttribute(dataset, "sigmay", sigmaydata, nt, ny);

    if(sxdata)
        createAndWriteDouble3DAttribute(dataset, "sx", sxdata, nt, nx, np);
    if(sydata)
        createAndWriteDouble3DAttribute(dataset, "sy", sydata, nt, ny, np);
    if(ssigmaydata)
        createAndWriteDouble3DAttribute(dataset, "ssigmay", ssigmaydata, nt, ny, np);
    // TODO: sssigmaz

    H5Fclose(file_id);

//    /** time derivative */
//    double *am_xdotdata;
//    /** parameter derivative of time derivative */
//    double *am_dxdotdpdata;
//    /** state derivative of observables */
//    double *am_dydxdata;
//    /** parameter derivative of observables */
//    double *am_dydpdata;
//    /** Jacobian of differential equation right hand side */
//    double *am_Jdata;
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
//    /** state */
//    double *am_xdata;
//    /** parameter derivative of state */
//    double *am_sxdata;
//    /** observable */
//    double *am_ydata;
//    /** observable standard deviation */
//    double *am_sigmaydata;
//    /** parameter derivative of observable */
//    double *am_sydata;
//    /** parameter derivative of observable standard deviation */
//    double *am_ssigmaydata;

//    /** number of integration steps forward problem */
//    double *am_numstepsdata;
//    /** number of integration steps backward problem */
//    double *am_numstepsSdata;
//    /** number of right hand side evaluations forward problem */
//    double *am_numrhsevalsdata;
//    /** number of right hand side evaluations backwad problem */
//    double *am_numrhsevalsSdata;
//    /** employed order forward problem */
//    double *am_orderdata;

//    /** likelihood value */
//    double *am_llhdata;
//    /** chi2 value */
//    double *am_chi2data;
//    /** parameter derivative of likelihood */
//    double *am_sllhdata;
//    /** second order parameter derivative of likelihood */
//    double *am_s2llhdata;

}


double getDoubleScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName) {
    double doubleScalar;
    H5LTget_attribute_double(file_id, optionsObject, attributeName, &doubleScalar);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %e\n", attributeName, doubleScalar);
#endif
    return doubleScalar;
}

int getIntScalarAttribute(hid_t file_id, const char* optionsObject, const char* attributeName) {
    int intScalar;
    H5LTget_attribute_int(file_id, optionsObject, attributeName, &intScalar);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d\n", attributeName, intScalar);
#endif
    return intScalar;
}

void getDoubleArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;
    herr_t status;
    status = H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);
    if(status < 0) {
        fprintf(stderr, "Error in getDoubleArrayAttribute: Cannot read attribute '%s' of '%s'\n", attributeName, optionsObject);
        #ifdef AMI_HDF5_H_DEBUG
        void *array[10];
        size_t size;
        size = backtrace(array, 10);
        backtrace_symbols_fd(array, size, STDERR_FILENO);
        #endif
    }
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d: ", attributeName, *length);
#endif
    *destination = new double[*length]; // vs. type_size
    status = H5LTget_attribute_double(file_id, optionsObject, attributeName, *destination);
    if(status < 0)
        fprintf(stderr, "Error in getDoubleArrayAttribute: Cannot read attribute '%s' of '%s'\n", attributeName, optionsObject);

#ifdef AMI_HDF5_H_DEBUG
    printfArray(*destination, *length, "%e ");
    printf("\n");
#endif
}

void getDoubleArrayAttribute2D(hid_t file_id, const char* optionsObject, const char* attributeName, double **destination, hsize_t *m, hsize_t *n) {
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
        getDoubleArrayAttribute(file_id, optionsObject, attributeName, destination, m);
    } else {
#ifdef AMI_HDF5_H_DEBUG
        printf("%s: %d x %d: ", attributeName, dims[0], dims[1]);
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

void getIntArrayAttribute(hid_t file_id, const char* optionsObject, const char* attributeName, int **destination, hsize_t *length) {
    H5T_class_t type_class;
    size_t type_size;

    H5LTget_attribute_info(file_id, optionsObject, attributeName, length, &type_class, &type_size);
#ifdef AMI_HDF5_H_DEBUG
    printf("%s: %d: ", attributeName, *length);
#endif
    *destination = new int[*length];
    H5LTget_attribute_int(file_id, optionsObject, attributeName, *destination);
}

// option for overwrite
void createAndWriteDouble2DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n) {
    const hsize_t adims[] = {m, n};

    if(H5Aexists(dataset, attributeName) > 0) {
        H5Adelete(dataset, attributeName);
    }

    hid_t space = H5Screate_simple(2, adims, NULL);
    hid_t attr = H5Acreate2(dataset, attributeName, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr, H5T_NATIVE_DOUBLE, buffer);
}

void createAndWriteDouble3DAttribute(hid_t dataset, const char *attributeName, const double *buffer, hsize_t m, hsize_t n, hsize_t o) {
    const hsize_t adims[] = {m, n, o};

    if(H5Aexists(dataset, attributeName) > 0) {
        H5Adelete(dataset, attributeName);
    }

    hid_t space = H5Screate_simple(3, adims, NULL);
    hid_t attr = H5Acreate2(dataset, attributeName, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr, H5T_NATIVE_DOUBLE, buffer);
}


void setAttributeIntFromDouble(hid_t file_id, const char *obj_name, const char *attr_name, const double *bufferDouble, size_t size)
{
    int intBuffer[size];
    for(int i = 0; i < size; ++i)
        intBuffer[i] = bufferDouble[i];

    H5LTset_attribute_int(file_id, obj_name, attr_name, intBuffer, size);
}
