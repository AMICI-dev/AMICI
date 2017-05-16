#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupJakstatAdjointO2)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(groupJakstatAdjointO2, testSensitivity2) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/sensi2forward/options");
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_jakstat_adjoint/sensi2forward/data");
    int status;
    CHECK_FALSE(edata == NULL);
    CHECK_FALSE(udata == NULL);

    ReturnData *rdata = getSimulationResults(udata, edata, &status);

    CHECK_EQUAL(0, status);

    verifyReturnData("/model_jakstat_adjoint/sensi2forward/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    freeReturnData(rdata);
    freeExpData(edata);
    freeUserData(udata);
}


