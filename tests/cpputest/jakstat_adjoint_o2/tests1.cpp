#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>

#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupJakstatAdjointO2)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(groupJakstatAdjointO2, testSensitivityForward2) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_jakstat_adjoint/sensi2forward/");
    delete model;
}

TEST(groupJakstatAdjointO2, testSensitivityAdjoint2) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_jakstat_adjoint/sensi2adjoint/");
    delete model;
}




