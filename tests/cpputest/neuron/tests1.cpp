#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupNeuron)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupNeuron, testSimulation) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_neuron/nosensi/options");
    ExpData *edata = NULL;

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_neuron/nosensi/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    delete rdata;
    delete udata;
}

TEST(groupNeuron, testSensitivityForward) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_neuron/sensiforward/options");
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_neuron/sensiforward/data");

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_neuron/sensiforward/results", rdata, udata, TEST_ATOL, TEST_RTOL);

    delete rdata;
    delete edata;
    delete udata;
}

