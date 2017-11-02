#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.h>
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>
#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupRobertson)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupRobertson, testSimulation) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_robertson/nosensi/");
    delete model;
}

TEST(groupRobertson, testSimulationExpData) {

}

TEST(groupRobertson, testSensitivityForward) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_robertson/sensiforward/");
    delete model;
}

TEST(groupRobertson, testSensitivityState) {

}

TEST(groupRobertson, testSensitivityAdjoint) {

}


