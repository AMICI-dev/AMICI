#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

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
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_steadystate/nosensi/options", model);
    ExpData *edata = NULL;

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_steadystate/nosensi/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);

    delete model;
    delete rdata;
    delete udata;
}

TEST(groupSteadystate, testSensitivityForward) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_steadystate/sensiforward/options", model);
    ExpData *edata = NULL;

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_steadystate/sensiforward/results", rdata, udata, model, TEST_ATOL, TEST_RTOL);

    delete model;
    delete rdata;
    delete udata;
}

