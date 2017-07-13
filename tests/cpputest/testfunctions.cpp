#include "testfunctions.h"
#include "include/ami_hdf5.h"
#include "include/amici_interface_cpp.h"
#include <cstring>
#include <execinfo.h>
#include <cstdio>
#include <unistd.h>

ExpData *getTestExpData(const UserData *udata) {
    ExpData *edata = new ExpData(udata);
    return edata;
}

bool withinTolerance(double expected, double actual, double atol, double rtol, int index) {
    bool withinTol =  fabs(expected - actual) <= atol || fabs((expected - actual) / (rtol + expected)) <= rtol;

    if(!withinTol) {
        fprintf(stderr, "ERROR: Expected value %e, but was %e at index %d.\n",expected, actual, index);
        fprintf(stderr, "       Relative error: %e (tolerance was %e)\n", fabs((expected - actual) / (rtol + expected)), rtol);
        fprintf(stderr, "       Absolute error: %e (tolerance was %e)\n", fabs(expected - actual), atol);
        printBacktrace(12);
    }

    return withinTol;
}

void checkEqualArray(const double *expected, const double *actual, int length, double atol, double rtol) {
    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i], actual[i], atol, rtol, i);
        CHECK_TRUE(withinTol);
    }
}

void checkEqualArrayStrided(const double *expected, const double *actual, int length, int strideExpected, int strideActual, double atol, double rtol) {
    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i * strideExpected], actual[i * strideActual], atol, rtol, i);
        CHECK_TRUE(withinTol);
    }
}

void verifyReturnData(const char* resultPath, const ReturnData *rdata, const UserData *udata, double atol, double rtol) {
    CHECK_FALSE(udata == NULL);
    CHECK_FALSE(rdata == NULL);

    // compare to saved data in hdf file
    hid_t file_id = H5Fopen(HDFFILE, H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t m, n, o;
    double *expected;

    double llhExp = AMI_HDF5_getDoubleScalarAttribute(file_id, resultPath, "llh");
    // TODO: ignores Inf and NaN results; need to check with format in HDF5
    if(! isinf(*rdata->llh) || isnan(*rdata->llh))
        CHECK_TRUE(withinTolerance(llhExp, *rdata->llh, atol, rtol, 1));

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "x", &expected, &m, &n);
    checkEqualArray(expected, rdata->x, udata->nt * udata->nxtrue, atol, rtol);
    delete[] expected;

//    CHECK_EQUAL(AMICI_O2MODE_FULL, udata->o2mode);

    if(AMI_HDF5_attributeExists(file_id, resultPath, "J")) {
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "J", &expected, &m, &n);
        checkEqualArray(expected, rdata->J, udata->nx * udata->nx, atol, rtol);
        delete[] expected;
    }

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "y", &expected, &m, &n);
    checkEqualArray(expected, rdata->y, udata->nt * udata->nytrue, atol, rtol);
    delete[] expected;
    
    if(udata->nz>0) {
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "z", &expected, &m, &n);
        checkEqualArray(expected, rdata->z, udata->nmaxevent * udata->nztrue, atol, rtol);
        delete[] expected;
        
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "rz", &expected, &m, &n);
        checkEqualArray(expected, rdata->rz, udata->nmaxevent * udata->nztrue, atol, rtol);
        delete[] expected;
        
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "sigmaz", &expected, &m, &n);
        checkEqualArray(expected, rdata->sigmaz, udata->nmaxevent * udata->nztrue, atol, rtol);
        delete[] expected;
    }


    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "xdot", &expected, &m, &n);
    checkEqualArray(expected, rdata->xdot, udata->nxtrue, atol, rtol);
    delete[] expected;

    if(udata->sensi >= AMICI_SENSI_ORDER_FIRST) {
        verifyReturnDataSensitivities(file_id, resultPath, rdata, udata, atol, rtol);
    } else {
        POINTERS_EQUAL(NULL, rdata->sllh);
        POINTERS_EQUAL(NULL, rdata->s2llh);
    }

    H5Fclose(file_id);
}

void verifyReturnDataSensitivities(hid_t file_id, const char* resultPath, const ReturnData *rdata, const UserData*udata, double atol, double rtol) {
    hsize_t m, n, o;
    double *expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "sllh", &expected, &m, &n);
    checkEqualArray(expected, rdata->sllh, udata->np, atol, rtol);
    delete[] expected;

    if(udata->sensi_meth == AMICI_SENSI_FSA) {
        AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "sx", &expected, &m, &n, &o);
        for(int ip = 0; ip < udata->nplist; ++ip)
            checkEqualArray(&expected[ip * udata->nt * udata->nxtrue],
                    &rdata->sx[ip * udata->nt * udata->nx],
                    udata->nt * udata->nxtrue, atol, rtol);
        delete[] expected;

        AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "sy", &expected, &m, &n, &o);
        for(int ip = 0; ip < udata->nplist; ++ip)
            checkEqualArray(&expected[ip * udata->nt * udata->nytrue],
                    &rdata->sy[ip * udata->nt * udata->ny],
                    udata->nt * udata->nytrue, atol, rtol);
        delete[] expected;
        
        if(udata->nz>0) {
            AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "sz", &expected, &m, &n, &o);
            for(int ip = 0; ip < udata->nplist; ++ip)
                checkEqualArray(&expected[ip * udata->nmaxevent * udata->nztrue],
                                &rdata->sz[ip * udata->nmaxevent * udata->nz],
                                udata->nmaxevent * udata->nztrue, atol, rtol);
            delete[] expected;
            
            AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "srz", &expected, &m, &n, &o);
            for(int ip = 0; ip < udata->nplist; ++ip)
                checkEqualArray(&expected[ip * udata->nmaxevent * udata->nztrue],
                                &rdata->srz[ip * udata->nmaxevent * udata->nz],
                                udata->nmaxevent * udata->nztrue, atol, rtol);
            delete[] expected;
        }
    }

    AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "ssigmay", &expected, &m, &n, &o);
    for(int ip = 0; ip < udata->nplist; ++ip)
        checkEqualArray(&expected[ip * udata->nt * udata->nytrue],
                &rdata->ssigmay[ip * udata->nt * udata->ny],
                udata->nt * udata->nytrue, atol, rtol);
    delete[] expected;
    
    if(udata->nz>0) {
        AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "ssigmaz", &expected, &m, &n, &o);
        for(int ip = 0; ip < udata->nplist; ++ip)
            checkEqualArray(&expected[ip * udata->nmaxevent * udata->nztrue],
                            &rdata->ssigmaz[ip * udata->nmaxevent * udata->nz],
                            udata->nmaxevent * udata->nztrue, atol, rtol);
        delete[] expected;
    }

    if(udata->sensi >= AMICI_SENSI_ORDER_SECOND) {
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "s2llh", &expected, &m, &n);
        checkEqualArray(expected, rdata->s2llh, (udata->nJ-1) * udata->nplist, atol, rtol);
        delete[] expected;
    } else {
        POINTERS_EQUAL(NULL, rdata->s2llh);
        POINTERS_EQUAL(NULL, rdata->s2rz);
    }

}


void printBacktrace(int depth) {
    void *array[depth];
    size_t size;
    size = backtrace(array, depth);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
}
