#include "testfunctions.h"

#include "wrapfunctions.h"
#include <cstring>

#include <gtest/gtest.h>

TEST(ExampleSteadystate, Default)
{
    amici::simulateWithDefaultOptions();
}

TEST(ExampleSteadystate, ModelFromHDF5)
{
    // Test reading some python-written options
    std::vector<double> pExp{ 1, 0.5, 0.4, 2, 0.1 };
    std::vector<double> kExp{ 0.1, 0.4, 0.7, 1.0 };

    auto model = amici::generic_model::getModel();
    amici::hdf5::readModelDataFromHDF5(
      NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");

    amici::checkEqualArray(
      kExp, model->getFixedParameters(), TEST_ATOL, TEST_RTOL, "k");
    ASSERT_EQ(51, model->nt());
    ASSERT_EQ(0.0, model->getTimepoint(0));
    ASSERT_EQ(100.0, model->getTimepoint(model->nt() - 2));
    ASSERT_EQ(INFINITY, model->getTimepoint(model->nt() - 1));

    for (int i = 0; i < model->np(); ++i) {
        ASSERT_EQ(pExp[i], model->getUnscaledParameters()[i]);
        ASSERT_EQ(log10(pExp[i]), model->getParameters()[i]);
    }
}

TEST(ExampleSteadystate, Inequality)
{
    auto modelA = amici::generic_model::getModel();
    auto modelB = std::make_unique<amici::Model_Test>();

    ASSERT_FALSE(*modelA == *modelB);
}

TEST(ExampleSteadystate, CopyModel)
{
    auto modelA = amici::generic_model::getModel();
    auto modelB = std::unique_ptr<amici::Model>(modelA->clone());

    ASSERT_EQ(*modelA, *modelB);
}

TEST(ExampleSteadystate, CloneModel)
{
    auto modelA = amici::generic_model::getModel();
    auto modelB = std::make_unique<
        amici::model_model_steadystate::Model_model_steadystate>();

    ASSERT_EQ(*modelA, *modelB);
}

TEST(ExampleSteadystate, ExpDataFromReturnData)
{
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();

    amici::hdf5::readModelDataFromHDF5(
      NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(
      NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    auto rdata = runAmiciSimulation(*solver, nullptr, *model);
    auto edata = amici::ExpData(*rdata, 0.1, 0.1);
    runAmiciSimulation(*solver, &edata, *model);
}

TEST(ExampleSteadystate, ReuseSolver)
{
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();

    amici::hdf5::readModelDataFromHDF5(
      NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(
      NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    runAmiciSimulation(*solver, nullptr, *model);
    runAmiciSimulation(*solver, nullptr, *model);
}

TEST(ExampleSteadystate, Rethrow)
{
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();

    amici::hdf5::readModelDataFromHDF5(
      NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(
      NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    // p = NaN will raise amici::IntegrationFailure
    auto p = model->getParameters();
    std::fill(p.begin(), p.end(), std::nan(""));
    model->setParameters(p);

    // must not throw
    runAmiciSimulation(*solver, nullptr, *model);

    // must throw
    ASSERT_THROW(runAmiciSimulation(*solver, nullptr, *model, true),
                 amici::IntegrationFailure);
}


TEST(ExampleSteadystate, Maxtime)
{
    auto model = amici::generic_model::getModel();
    auto solver = model->getSolver();

    amici::hdf5::readModelDataFromHDF5(
        NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::readSolverSettingsFromHDF5(
        NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    auto rdata = runAmiciSimulation(*solver, nullptr, *model);
    ASSERT_EQ(amici::AMICI_SUCCESS, rdata->status);

    solver->setMaxTime(0.000001);

    // must throw
    rdata = runAmiciSimulation(*solver, nullptr, *model);
    ASSERT_EQ(amici::AMICI_MAX_TIME_EXCEEDED, rdata->status);
}

TEST(ExampleSteadystate, InitialStatesNonEmpty)
{
    auto model = amici::generic_model::getModel();
    ASSERT_FALSE(model->getInitialStates().empty());
}

TEST(ExampleSteadystate, InitialStateSensitivitiesNonEmpty)
{
    auto model = amici::generic_model::getModel();
    ASSERT_FALSE(model->getInitialStateSensitivities().empty());
}

TEST(ExampleSteadystate, Simulation)
{
    amici::simulateVerifyWrite(
      "/model_steadystate/nosensi/", 100 * TEST_ATOL, 100 * TEST_RTOL);
}

TEST(ExampleSteadystate, SensitivityForward)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiforward/");
}

TEST(ExampleSteadystate, SensitivityForwardPlist)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiforwardplist/");
}

TEST(ExampleSteadystate, SensitivityForwardErrorInt)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiforwarderrorint/");
}

TEST(ExampleSteadystate, SensitivityForwardErrorNewt)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiforwarderrornewt/");
}

TEST(ExampleSteadystate, SensitivityForwardDense)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiforwarddense/");
}

TEST(ExampleSteadystate, SensitivityForwardSPBCG)
{
    amici::simulateVerifyWrite(
      "/model_steadystate/nosensiSPBCG/", 10 * TEST_ATOL, 10 * TEST_RTOL);
}

TEST(ExampleSteadystate, SensiFwdNewtonPreeq)
{
    amici::simulateVerifyWrite("/model_steadystate/sensifwdnewtonpreeq/");
}

TEST(ExampleSteadystate, SensiAdjNewtonPreeq)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiadjnewtonpreeq/");
}

TEST(ExampleSteadystate, SensiFwdSimPreeq)
{
    amici::simulateVerifyWrite("/model_steadystate/sensifwdsimpreeq/");
}

TEST(ExampleSteadystate, SensiFwdSimPreeqFSA)
{
    amici::simulateVerifyWrite("/model_steadystate/sensifwdsimpreeqFSA/");
}

TEST(ExampleSteadystate, SensiAdjSimPreeq)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiadjsimpreeq/");
}

TEST(ExampleSteadystate, SensiAdjSimPreeqFSA)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiadjsimpreeqFSA/");
}

TEST(ExampleSteadystate, SensiFwdByhandPreeq)
{
    amici::simulateVerifyWrite("/model_steadystate/sensifwdbyhandpreeq/");
}

TEST(ExampleSteadystate, SensiAdjByhandPreeq)
{
    amici::simulateVerifyWrite("/model_steadystate/sensiadjbyhandpreeq/");
}
