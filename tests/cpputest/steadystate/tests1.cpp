#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupSteadystate)
{
    void setup() {

    }

    void teardown() {

    }
};



TEST(groupSteadystate, testSimulation) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_steadystate/nosensi/options");
    ExpData *edata = NULL;

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    verifyReturnData("/model_steadystate/nosensi/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    freeReturnData(rdata);
    freeExpData(edata);
    delete udata;
}

TEST(groupSteadystate, testSensitivityForward) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_steadystate/sensiforward/options");
    ExpData *edata = NULL;

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    verifyReturnData("/model_steadystate/sensiforward/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    freeReturnData(rdata);
    freeExpData(edata);
    delete udata;
}

