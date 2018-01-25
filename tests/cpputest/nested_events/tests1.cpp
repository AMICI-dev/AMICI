#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupEvents)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupEvents, testSimulation) {
    amici::simulateAndVerifyFromFile("/model_nested_events/nosensi/");
}

TEST(groupEvents, testSensitivityForward) {
    amici::simulateAndVerifyFromFile("/model_nested_events/sensiforward/");
}


