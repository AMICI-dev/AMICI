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
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_events/nosensi/");
    delete model;
}

TEST(groupEvents, testSensitivityForward) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_events/sensiforward/");
    delete model;
}


