#include "testfunctions.h"
#include <cstring>
ExpData *getTestExpData() {
    ExpData *edata = new ExpData;
    memset(edata, 0, sizeof(*edata));

    return edata;
}

bool withinTolerance(double expected, double actual, double atol, double rtol) {
    return fabs(expected - actual) <= atol || fabs((expected - actual) / expected) <= rtol;
}

void checkEqualArray(const double *expected, const double *actual, int length, double atol, double rtol) {
    for(int i = 0; i < length; ++i)
    {
        bool withinTol = withinTolerance(expected[i], actual[i], atol, rtol);

#ifndef __APPLE__
        if(!withinTol)
            std::cout<<i<<"/"<<length<<" "<<expected[i]<<" "<<actual[i]<<std::endl;
#endif

        CHECK_TRUE(withinTol);
    }
}

void verifyReturnData(const char* resultPath, const ReturnData *rdata, const UserData*udata, double atol, double rtol) {
    // compare to saved data in hdf file
    hid_t file_id = H5Fopen(HDFFILE, H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t m, n, o;

    double *expected;

    double llhExp = AMI_HDF5_getDoubleScalarAttribute(file_id, resultPath, "llh");
    // TODO: ignores Inf and NaN results; need to check with format in HDF5
    if(! isinf(*rdata->am_llhdata) || isnan(*rdata->am_llhdata))
        CHECK_TRUE(withinTolerance(llhExp, *rdata->am_llhdata, atol, rtol));

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "x", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_xdata, udata->am_nt * udata->am_nx, atol, rtol);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "J", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_Jdata, udata->am_nx * udata->am_nx, atol, rtol);
    delete[] expected;

//    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "numrhsevals", &expected, &m, &n);
//    checkEqualArray(expected, rdata->am_numrhsevalsdata, udata->am_nt, epsilon, blab);
//    delete[] expected;

//    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "numsteps", &expected, &m, &n);
//    checkEqualArray(expected, rdata->am_numstepsdata, udata->am_nt, epsilon, blab);
//    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "order", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_orderdata, udata->am_nt, atol, rtol);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "y", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_ydata, udata->am_nt, atol, rtol);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "sigmay", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_sigmaydata, udata->am_nt, atol, rtol);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "xdot", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_xdotdata, udata->am_nx, atol, rtol);
    delete[] expected;

    if(udata->am_sensi >= AMI_SENSI_ORDER_FIRST) {
        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "sllh", &expected, &m, &n);
        checkEqualArray(expected, rdata->am_sllhdata, udata->am_np, atol, rtol);
        delete[] expected;

//        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "numrhsevalsS", &expected, &m, &n);
//        checkEqualArray(expected, rdata->am_numrhsevalsSdata, udata->am_nt, epsilon, blab);
//        delete[] expected;

//        AMI_HDF5_getDoubleArrayAttribute2D(file_id, resultPath, "numstepsS", &expected, &m, &n);
//        checkEqualArray(expected, rdata->am_numstepsSdata, udata->am_nt, epsilon, blab);
//        delete[] expected;

        AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "ssigmay", &expected, &m, &n, &o);
        checkEqualArray(expected, rdata->am_ssigmaydata, udata->am_nt * udata->am_ny * udata->am_nplist, atol, rtol);
        delete[] expected;

        if(udata->am_sensi_meth == AMI_SENSI_FSA) {
            AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "sx", &expected, &m, &n, &o);
            checkEqualArray(expected, rdata->am_sxdata, udata->am_nt * udata->am_nx * udata->am_nplist, atol, rtol);
            delete[] expected;

            AMI_HDF5_getDoubleArrayAttribute3D(file_id, resultPath, "sy", &expected, &m, &n, &o);
            checkEqualArray(expected, rdata->am_sydata, udata->am_nt * udata->am_ny * udata->am_nplist, atol, rtol);
            delete[] expected;

        }

    } else {
        POINTERS_EQUAL(NULL, rdata->am_sllhdata);
        POINTERS_EQUAL(NULL, rdata->am_numrhsevalsSdata);
        POINTERS_EQUAL(NULL, rdata->am_numstepsSdata);

    }

    H5Fclose(file_id);
}
