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
    ExpData *edata = getTestExpData();

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    // TODO proper paths /teststeadystate1/...
    verifyReturnData("/model_steadystate/nosensi/results", rdata, udata);

    freeReturnData(rdata);
    freeExpData(edata);
    freeUserData(udata);
}



