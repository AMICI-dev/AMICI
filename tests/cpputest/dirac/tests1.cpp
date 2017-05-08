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


/*
 * Test for mem leaks in UserData initalization / destruction
 */
TEST(groupDirac, testCreateAndFreeUserData) {
    UserData *udata = getTestUserData();

    freeUserData(udata);
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
    verifyReturnData("/model_dirac/nosensi/results", rdata, udata);

    freeReturnData(rdata);
    freeExpData(edata);
    freeUserData(udata);
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
    verifyReturnData("/model_dirac/sensiforward/results", rdata, udata);

    freeReturnData(rdata);
    freeExpData(edata);
    freeUserData(udata);
}

TEST(groupDirac, testSensitivityState) {

}

TEST(groupDirac, testSensitivityAdjoint) {

}


