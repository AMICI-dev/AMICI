#include "testfunctions.h"

#include "wrapfunctions.h"
#include "model_steadystate_py.h"

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

    auto model = amici::generic_model::get_model();
    amici::hdf5::read_model_data_from_hdf5(
      NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");

    amici::checkEqualArray(
      kExp, model->get_fixed_parameters(), TEST_ATOL, TEST_RTOL, "k");
    ASSERT_EQ(51, model->nt());
    ASSERT_EQ(0.0, model->get_timepoint(0));
    ASSERT_EQ(100.0, model->get_timepoint(model->nt() - 2));
    ASSERT_EQ(INFINITY, model->get_timepoint(model->nt() - 1));

    for (int i = 0; i < model->np(); ++i) {
        ASSERT_EQ(pExp[i], model->get_unscaled_parameters()[i]);
        ASSERT_EQ(log10(pExp[i]), model->get_free_parameters()[i]);
    }
}

TEST(ExampleSteadystate, Inequality)
{
    auto modelA = amici::generic_model::get_model();
    auto modelB = std::make_unique<amici::Model_Test>();

    ASSERT_FALSE(*modelA == *modelB);
}

TEST(ExampleSteadystate, CopyModel)
{
    auto modelA = amici::generic_model::get_model();
    auto modelB = std::unique_ptr<amici::Model>(modelA->clone());

    ASSERT_EQ(*modelA, *modelB);
}

TEST(ExampleSteadystate, CloneModel)
{
    auto modelA = amici::generic_model::get_model();
    auto modelB = std::make_unique<
        amici::model_model_steadystate_py::Model_model_steadystate_py>();

    ASSERT_EQ(*modelA, *modelB);
}

TEST(ExampleSteadystate, ExpDataFromReturnData)
{
    auto model = amici::generic_model::get_model();
    auto solver = model->create_solver();

    amici::hdf5::read_model_data_from_hdf5(
      NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::read_solver_settings_from_hdf5(
      NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    auto rdata = run_simulation(*solver, nullptr, *model);
    auto edata = amici::ExpData(*rdata, 0.1, 0.1);
    run_simulation(*solver, &edata, *model);
}

TEST(ExampleSteadystate, ReuseSolver)
{
    auto model = amici::generic_model::get_model();
    auto solver = model->create_solver();

    amici::hdf5::read_model_data_from_hdf5(
      NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::read_solver_settings_from_hdf5(
      NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    run_simulation(*solver, nullptr, *model);
    run_simulation(*solver, nullptr, *model);
}

TEST(ExampleSteadystate, Rethrow)
{
    auto model = amici::generic_model::get_model();
    auto solver = model->create_solver();

    amici::hdf5::read_model_data_from_hdf5(
      NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::read_solver_settings_from_hdf5(
      NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    // p = NaN will raise amici::IntegrationFailure
    auto p = model->get_free_parameters();
    std::fill(p.begin(), p.end(), std::nan(""));
    model->set_free_parameters(p);

    // must not throw
    run_simulation(*solver, nullptr, *model);

    // must throw
    ASSERT_THROW(run_simulation(*solver, nullptr, *model, true),
                 amici::IntegrationFailure);
}


TEST(ExampleSteadystate, Maxtime)
{
    auto model = amici::generic_model::get_model();
    auto solver = model->create_solver();

    amici::hdf5::read_model_data_from_hdf5(
        NEW_OPTION_FILE, *model, "/model_steadystate/nosensi/options");
    amici::hdf5::read_solver_settings_from_hdf5(
        NEW_OPTION_FILE, *solver, "/model_steadystate/nosensi/options");

    // Ensure the solver needs sufficiently many steps that time is checked
    // at least once during integration
    solver->set_relative_tolerance(1e-14);

    auto rdata = run_simulation(*solver, nullptr, *model);
    ASSERT_EQ(amici::AMICI_SUCCESS, rdata->status);

    solver->set_max_time(0.000001);

    // must throw
    rdata = run_simulation(*solver, nullptr, *model);
    ASSERT_EQ(amici::AMICI_MAX_TIME_EXCEEDED, rdata->status);
}

TEST(ExampleSteadystate, InitialStatesNonEmpty)
{
    auto model = amici::generic_model::get_model();
    ASSERT_FALSE(model->get_initial_state().empty());
}

TEST(ExampleSteadystate, InitialStateSensitivitiesNonEmpty)
{
    auto model = amici::generic_model::get_model();
    ASSERT_FALSE(model->get_initial_state_sensitivities().empty());
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
    amici::simulateVerifyWrite("/model_steadystate/sensiforwarddense/",
                               1e-9, TEST_RTOL);
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
