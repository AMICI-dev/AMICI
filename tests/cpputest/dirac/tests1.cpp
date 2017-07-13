#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.h>
#include <include/ami_hdf5.h>
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
    UserData *udata = new UserData(getUserData());

    udata->qpositivex = new double[udata->nx];
    udata->p = new double[udata->np];
    udata->k = new double[udata->nk];
    udata->idlist = new realtype[udata->np]();
    for(int i = 0; i < udata->np; ++i)
        udata->idlist[i] = 0;

    udata->initTemporaryFields();

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
    UserData udata(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, AMICI_SCALING_NONE, AMICI_O2MODE_NONE);
    
    ExpData *edata = getTestExpData(&udata);
    
    delete edata;
}

/*
 * Test for mem leaks in ReturnData initalization / destruction
 */

TEST(groupDirac, testCreateAndFreeReturnData) {
    UserData udata(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, AMICI_SCALING_NONE, AMICI_O2MODE_NONE);
    ReturnData *rdata = new ReturnData(&udata);

    delete rdata;
}


TEST(groupDirac, testSimulation) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_dirac/nosensi/options");
    ExpData *edata = NULL;

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    // TODO proper paths /testDirac1/...
    verifyReturnData("/model_dirac/nosensi/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    delete rdata;
    delete udata;
}

TEST(groupDirac, testSimulationExpData) {

}

TEST(groupDirac, testSensitivityForward) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_dirac/sensiforward/options");
    ExpData *edata = NULL;

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    // TODO proper paths /testDirac1/...
    verifyReturnData("/model_dirac/sensiforward/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    delete rdata;
    delete udata;
}

TEST(groupDirac, testSensitivityState) {

}

TEST(groupDirac, testSensitivityAdjoint) {

}


