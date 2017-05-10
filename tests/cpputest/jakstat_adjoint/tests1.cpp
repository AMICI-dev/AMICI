#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <cstring>
#include "wrapfunctions.h"

#define TEST_EPSILON 1e-10

TEST_GROUP(groupJakstatAdjoint)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(groupJakstatAdjoint, testSimulation) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/nosensi/options");
    ExpData *edata = NULL;

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    verifyReturnData("/model_jakstat_adjoint/nosensi/results", rdata, udata, TEST_EPSILON);

    freeReturnData(rdata);
    freeExpData(edata);
    freeUserData(udata);
}

TEST(groupJakstatAdjoint, testSensitivityForward) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/sensiforward/options");
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_jakstat_adjoint/sensiforward/data");

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    verifyReturnData("/model_jakstat_adjoint/sensiforward/results", rdata, udata, TEST_EPSILON);

    freeReturnData(rdata);
    freeExpData(edata);
    freeUserData(udata);
}

TEST(groupJakstatAdjoint, testSensitivityAdjoint) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_jakstat_adjoint/sensiadjoint/options");
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_jakstat_adjoint/sensiadjoint/data");

    int status;
    ReturnData *rdata = getSimulationResults(udata, edata, &status);
    CHECK_EQUAL(0, status);

    verifyReturnData("/model_jakstat_adjoint/sensiadjoint/results", rdata, udata, TEST_EPSILON);

    freeReturnData(rdata);
    freeExpData(edata);
    freeUserData(udata);
}


