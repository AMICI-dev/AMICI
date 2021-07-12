#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

TEST(ExampleJakstatAdjoint, Simulation)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/nosensi/");
}

TEST(ExampleJakstatAdjoint, SensitivityForward)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiforward/");
}

TEST(ExampleJakstatAdjoint, SensitivityForwardLogParam)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiforwardlogparam/");
}

TEST(ExampleJakstatAdjoint, SensitivityAdjoint)
{
    amici::simulateVerifyWrite("/model_jakstat_adjoint/sensiadjoint/");
}

TEST(ExampleJakstatAdjoint, SensitivityForwardEmptySensInd)
{
    amici::simulateVerifyWrite(
      "/model_jakstat_adjoint/sensiforwardemptysensind/");
}

TEST(ExampleJakstatAdjoint, SensitivityAdjointEmptySensInd)
{
    amici::simulateVerifyWrite(
      "/model_jakstat_adjoint/sensiadjointemptysensind/");
}

TEST(ExampleJakstatAdjoint, DISABLED_SensitivityAdjointUnusedNanOutputs)
{
    /* UN-IGNORE ONCE THIS MODEL HAS BEEN IMPORTED VIA PYTHON INTERFACE */
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();
    amici::hdf5::readModelDataFromHDF5(
      NEW_OPTION_FILE, *model, "/model_jakstat_adjoint/sensiadjoint/options");
    amici::hdf5::readSolverSettingsFromHDF5(
      NEW_OPTION_FILE, *solver, "/model_jakstat_adjoint/sensiadjoint/options");
    auto edata = amici::hdf5::readSimulationExpData(
      NEW_OPTION_FILE, "/model_jakstat_adjoint/sensiadjoint/data", *model);

    // Set output parameter p[10] to NaN and remove respective measurements
    // -> gradient should still be finite

    auto p = model->getParameters();
    p[10] = NAN;
    model->setParameters(p);

    auto d = edata->getObservedData();
    for (int it = 0; it < edata->nt(); ++it) {
        for (int iy = 0; iy < edata->nytrue(); ++iy) {
            if (iy == 1)
                d[it * edata->nytrue() + iy] = NAN;
        }
    }
    edata->setObservedData(d);

    auto rdata = runAmiciSimulation(*solver, edata.get(), *model);

    for (int i = 0; i < model->nplist(); ++i)
        ASSERT_FALSE(std::isnan(rdata->sllh[i]));
}

TEST(ExampleJakstatAdjoint, SensitivityReplicates)
{
    // Check that we can handle replicates correctly

    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();
    amici::hdf5::readModelDataFromHDF5(
      NEW_OPTION_FILE, *model, "/model_jakstat_adjoint/sensiadjoint/options");
    amici::hdf5::readSolverSettingsFromHDF5(
      NEW_OPTION_FILE, *solver, "/model_jakstat_adjoint/sensiadjoint/options");
    amici::ExpData edata(*model);

    // No replicate, no sensi
    edata.setTimepoints({ 10.0 });
    auto d = edata.getObservedData();
    for (int it = 0; it < edata.nt(); ++it) {
        for (int iy = 0; iy < edata.nytrue(); ++iy) {
            if (iy == 0) {
                d[it * edata.nytrue() + iy] = 1.0;
            } else {
                d[it * edata.nytrue() + iy] = NAN;
            }
        }
    }
    edata.setObservedData(d);
    edata.setObservedDataStdDev(1.0);

    solver->setSensitivityOrder(amici::SensitivityOrder::none);
    auto rdata1 = runAmiciSimulation(*solver, &edata, *model);
    auto llh1 = rdata1->llh;

    // forward + replicates
    edata.setTimepoints({ 10.0, 10.0 });
    d = edata.getObservedData();
    for (int it = 0; it < edata.nt(); ++it) {
        for (int iy = 0; iy < edata.nytrue(); ++iy) {
            if (iy == 0) {
                d[it * edata.nytrue() + iy] = 1.0;
            } else {
                d[it * edata.nytrue() + iy] = NAN;
            }
        }
    }
    edata.setObservedData(d);
    edata.setObservedDataStdDev(1.0);

    solver->setSensitivityOrder(amici::SensitivityOrder::first);
    solver->setSensitivityMethod(amici::SensitivityMethod::forward);
    auto rdata2 = runAmiciSimulation(*solver, &edata, *model);
    auto llh2 = rdata2->llh;
    ASSERT_NEAR(2.0 * llh1, llh2, 1e-6);

    // adjoint + replicates
    solver->setSensitivityMethod(amici::SensitivityMethod::adjoint);
    auto rdata3 = runAmiciSimulation(*solver, &edata, *model);
    auto llh3 = rdata3->llh;
    ASSERT_NEAR(llh2, llh3, 1e-6);
}
