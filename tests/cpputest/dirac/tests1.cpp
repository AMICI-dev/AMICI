#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.h>
#include <include/amici_interface_cpp.h>
#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupDirac)
{
    void setup() {

    }

    void teardown() {

    }
};


UserData *getTestUserData() {
    UserData *udata = getUserData();

    udata->qpositivex = new double[udata->nx];
    udata->p = new double[udata->np];
    udata->k = new double[udata->nk];
    udata->idlist = new realtype[udata->np]();
    for(int i = 0; i < udata->np; ++i)
        udata->idlist[i] = 0;

    processUserData(udata);

    return udata;
}


/*
 * Test for mem leaks in UserData initalization / destruction
 */
TEST(groupDirac, testCreateAndFreeUserData) {
    UserData *udata = getTestUserData();

    delete udata;
}

/*
 * Test for mem leaks in ExpData initalization / destruction
 */

TEST(groupDirac, testCreateAndFreeExpData) {
    ExpData *edata = getTestExpData();

    freeExpData(edata);
}

/*
 * Test for mem leaks in ReturnData initalization / destruction
 */

TEST(groupDirac, testCreateAndFreeReturnData) {
    ReturnData *rdata = new ReturnData;
    memset(rdata, 0, sizeof(*rdata));

    freeReturnData(rdata);
}


TEST(groupDirac, testSimulation) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_dirac/nosensi/options");
    ExpData *edata = NULL;

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    // TODO proper paths /testDirac1/...
    verifyReturnData("/model_dirac/nosensi/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    freeReturnData(rdata);
    freeExpData(edata);
    delete udata;
}

TEST(groupDirac, testSimulationExpData) {

}

TEST(groupDirac, testSensitivityForward) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_dirac/sensiforward/options");
    ExpData *edata = NULL;

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    // TODO proper paths /testDirac1/...
    verifyReturnData("/model_dirac/sensiforward/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    freeReturnData(rdata);
    freeExpData(edata);
    delete udata;
}

TEST(groupDirac, testSensitivityState) {

}

TEST(groupDirac, testSensitivityAdjoint) {

}


