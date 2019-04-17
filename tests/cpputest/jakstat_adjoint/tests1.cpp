#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

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
    amici::simulateVerifyWrite("/model_jakstat_adjoint/nosensi/");
}

TEST(groupJakstatAdjoint, testSensitivityForward) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiforward/");
}

TEST(groupJakstatAdjoint, testSensitivityForwardLogParam) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiforwardlogparam/");
}

TEST(groupJakstatAdjoint, testSensitivityAdjoint) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiadjoint/");
}

TEST(groupJakstatAdjoint, testSensitivityForwardEmptySensInd) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiforwardemptysensind/");
}

TEST(groupJakstatAdjoint, testSensitivityAdjointEmptySensInd) {
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiadjointemptysensind/");
}

IGNORE_TEST(groupJakstatAdjoint, testSensitivityAdjointUnusedNanOutputs) {
    /* UN-IGNORE ONCE THIS MODEL HAS BEEN IMPORTED VIA PYTHON INTERFACE */
    auto model = getModel();
    auto solver = model->getSolver();
    amici::hdf5::readModelDataFromHDF5(
                NEW_OPTION_FILE, *model,
                "/model_jakstat_adjoint/sensiadjoint/options");
    amici::hdf5::readSolverSettingsFromHDF5(
                NEW_OPTION_FILE, *solver,
                "/model_jakstat_adjoint/sensiadjoint/options");
    auto edata = amici::hdf5::readSimulationExpData(
                NEW_OPTION_FILE, "/model_jakstat_adjoint/sensiadjoint/data",
                *model);

    // Set output parameter p[10] to NaN and remove respective measurements
    // -> gradient should still be finite

    auto p = model->getParameters();
    p[10] = NAN;
    model->setParameters(p);

    auto d = edata->getObservedData();
    for (int it = 0; it < edata->nt(); ++it) {
        for (int iy = 0; iy < edata->nytrue(); ++iy) {
            if(iy == 1)
                d[it * edata->nytrue() + it] = NAN;
        }
    }
    edata->setObservedData(d);

    auto rdata = runAmiciSimulation(*solver, edata.get(), *model);

    for(int i = 0; i < model->nplist(); ++i)
        CHECK_FALSE(std::isnan(rdata->sllh[i]));
}
