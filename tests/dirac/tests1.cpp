#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"
#include <iostream>
#include <include/amici.h>
#include <cstring>
#include "wrapfunctions.h"
#include "src/ami_hdf5.h"

#define HDFFILE "../expectedResults.h5"
#define TEST_EPSILON 1e-12
TEST_GROUP(group1)
{
    void setup() {

    }

    void teardown() {

    }
};


UserData *getTestUserData() {
    UserData *udata = new UserData;
    memset(udata, 0, sizeof(*udata));

    init_modeldims(udata);

    udata->am_qpositivex = new double[udata->am_nx];
    udata->am_p = new double[udata->am_np];
    udata->am_k = new double[udata->am_nk];
    udata->am_idlist = new realtype[udata->am_np]();
    for(int i = 0; i < udata->am_np; ++i)
        udata->am_idlist[i] = 0;

    processUserData(udata);

    return udata;
}

ExpData *getTestExpData() {
    ExpData *edata = new ExpData;
    memset(edata, 0, sizeof(*edata));

    return edata;
}


/*
 * Test for mem leaks in UserData initalization / destruction
 */
TEST(group1, testCreateAndFreeUserData) {
    UserData *udata = getTestUserData();

    freeUserData(udata);
}

/*
 * Test for mem leaks in ExpData initalization / destruction
 */

TEST(group1, testCreateAndFreeExpData) {
    ExpData *edata = getTestExpData();

    freeExpData(edata);
}

/*
 * Test for mem leaks in ReturnData initalization / destruction
 */

TEST(group1, testCreateAndFreeReturnData) {
    ReturnData *rdata = new ReturnData;
    memset(rdata, 0, sizeof(*rdata));

    freeReturnData(rdata);
}

void checkEqualArray(const double *expected, const double *actual, int length) {
    for(int i = 0; i < length; ++i) {
        // std::cout<<i<<"/"<<length<<" "<<expected[i]<<" "<<actual[i]<<std::endl;
        DOUBLES_EQUAL(expected[i], actual[i], TEST_EPSILON);
    }
}

void verifyReturnData(const char* hdfGroup, const ReturnData *rdata, const UserData*udata) {
    CHECK_TRUE(isinf(*rdata->am_llhdata) || isnan(*rdata->am_llhdata));

    // compare to saved data in hdf file
    hid_t file_id = H5Fopen(HDFFILE, H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t m, n;

    double *expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, "/results", "x", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_xdata, udata->am_nt * udata->am_nx);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, "/results", "J", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_Jdata, udata->am_nx * udata->am_nx);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, "/results", "numrhsevals", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_numrhsevalsdata, udata->am_nt);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, "/results", "numsteps", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_numstepsdata, udata->am_nt);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, "/results", "order", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_orderdata, udata->am_nt);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, "/results", "y", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_sigmaydata, udata->am_nt);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, "/results", "sigmay", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_ydata, udata->am_nt);
    delete[] expected;

    AMI_HDF5_getDoubleArrayAttribute2D(file_id, "/results", "xdot", &expected, &m, &n);
    checkEqualArray(expected, rdata->am_xdotdata, udata->am_nx);
    delete[] expected;

    H5Fclose(file_id);
}

TEST(group1, testSimulation) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/options");
    ExpData *edata = getTestExpData();

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    // TODO proper paths /testDirac1/...
    verifyReturnData("/", rdata, udata);

    freeReturnData(rdata);
    freeExpData(edata);
    freeUserData(udata);
}

TEST(group1, testSimulationExpData) {

}

TEST(group1, testSensitivityForward) {

}

TEST(group1, testSensitivityState) {

}

TEST(group1, testSensitivityAdjoint) {

}


