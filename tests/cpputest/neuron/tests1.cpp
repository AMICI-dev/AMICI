#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

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
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_neuron/nosensi/options", model);
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_neuron/sensiforward/data", model);

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_neuron/nosensi/results", rdata, udata, model, 10*TEST_ATOL, 10*TEST_RTOL);

    delete model;
    delete rdata;
    delete edata;
    delete udata;
}

TEST(groupNeuron, testSensitivityForward) {
    Model *model = getModel();

    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_neuron/sensiforward/options", model);
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_neuron/sensiforward/data", model);

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_neuron/sensiforward/results", rdata, udata, model, 10*TEST_ATOL, 10*TEST_RTOL);

    delete model;
    delete rdata;
    delete edata;
    delete udata;
}

