#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupJakstatAdjoint)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(groupJakstatAdjoint, testSimulation) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_jakstat_adjoint/nosensi/");
    delete model;
}

TEST(groupJakstatAdjoint, testSensitivityForward) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_jakstat_adjoint/sensiforward/");
    delete model;
}

TEST(groupJakstatAdjoint, testSensitivityForwardLogParam) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_jakstat_adjoint/sensiforwardlogparam/");
    delete model;
}

TEST(groupJakstatAdjoint, testSensitivityAdjoint) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_jakstat_adjoint/sensiadjoint/");
    delete model;
}

TEST(groupJakstatAdjointO2, testSensitivityAdjoint2) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_jakstat_adjoint/sensi2adjoint/");
    delete model;
}

