#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupNeuronO2)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(groupNeuronO2, testSensitivity2) {
    // read simulation options
    UserData *udata = AMI_HDF5_readSimulationUserDataFromFileName(HDFFILE, "/model_neuron/sensi2forward/options");
    ExpData *edata = AMI_HDF5_readSimulationExpData(HDFFILE, udata, "/model_neuron/sensi2forward/data");
    CHECK_FALSE(edata == NULL);
    CHECK_FALSE(udata == NULL);

    ReturnData *rdata = getSimulationResults(udata, edata);
    CHECK_EQUAL(0, *rdata->status);

    verifyReturnData("/model_neuron/sensi2forward/results", rdata, udata, 10*TEST_ATOL, 10*TEST_RTOL);

    delete rdata;
    delete edata;
    delete udata;
}



