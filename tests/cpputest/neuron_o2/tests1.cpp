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
    amici::simulateAndVerifyFromFile("/model_neuron/sensi2forward/", 10*TEST_ATOL, 10*TEST_RTOL);
}



