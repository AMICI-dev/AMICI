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

TEST(ExampleJakstatAdjoint, SensitivityAdjointUnusedNanOutputs)
{
    /* UN-IGNORE ONCE THIS MODEL HAS BEEN IMPORTED VIA PYTHON INTERFACE */
    auto model = amici::generic_model::get_model();
    auto solver = model->create_solver();
    amici::hdf5::read_model_data_from_hdf5(
      NEW_OPTION_FILE, *model, "/model_jakstat_adjoint/sensiadjoint/options");
    amici::hdf5::read_solver_settings_from_hdf5(
      NEW_OPTION_FILE, *solver, "/model_jakstat_adjoint/sensiadjoint/options");
    auto edata = amici::hdf5::read_exp_data_from_hdf5(
      NEW_OPTION_FILE, "/model_jakstat_adjoint/sensiadjoint/data", *model);

    // Set output parameter p[10] to NaN and remove respective measurements
    // -> gradient should still be finite

    auto p = model->get_free_parameters();
    p[10] = NAN;
    model->set_free_parameter_by_id("offset_tSTAT", NAN);

    auto iy = 1;
    Expects(model->get_observable_ids()[iy] == "obs_tSTAT");
    auto d = edata->get_measurements();
    for (int it = 0; it < edata->nt(); ++it) {
        d[it * edata->nytrue() + iy] = NAN;
    }
    edata->set_measurements(d);

    auto rdata = run_simulation(*solver, edata.get(), *model);

    for (int i = 0; i < model->nplist(); ++i)
        ASSERT_FALSE(std::isnan(rdata->sllh[i]));
}

TEST(ExampleJakstatAdjoint, SensitivityReplicates)
{
    // Check that we can handle replicates correctly

    auto model = amici::generic_model::get_model();
    auto solver = model->create_solver();
    amici::hdf5::read_model_data_from_hdf5(
      NEW_OPTION_FILE, *model, "/model_jakstat_adjoint/sensiadjoint/options");
    amici::hdf5::read_solver_settings_from_hdf5(
      NEW_OPTION_FILE, *solver, "/model_jakstat_adjoint/sensiadjoint/options");
    amici::ExpData edata(*model);

    // No replicate, no sensi
    edata.set_timepoints({ 10.0 });
    auto d = edata.get_measurements();
    for (int it = 0; it < edata.nt(); ++it) {
        for (int iy = 0; iy < edata.nytrue(); ++iy) {
            if (iy == 0) {
                d[it * edata.nytrue() + iy] = 1.0;
            } else {
                d[it * edata.nytrue() + iy] = NAN;
            }
        }
    }
    edata.set_measurements(d);
    edata.set_noise_scales(1.0);

    solver->set_sensitivity_order(amici::SensitivityOrder::none);
    auto rdata1 = run_simulation(*solver, &edata, *model);
    auto llh1 = rdata1->llh;

    // forward + replicates
    edata.set_timepoints({ 10.0, 10.0 });
    d = edata.get_measurements();
    for (int it = 0; it < edata.nt(); ++it) {
        for (int iy = 0; iy < edata.nytrue(); ++iy) {
            if (iy == 0) {
                d[it * edata.nytrue() + iy] = 1.0;
            } else {
                d[it * edata.nytrue() + iy] = NAN;
            }
        }
    }
    edata.set_measurements(d);
    edata.set_noise_scales(1.0);

    solver->set_sensitivity_order(amici::SensitivityOrder::first);
    solver->set_sensitivity_method(amici::SensitivityMethod::forward);
    auto rdata2 = run_simulation(*solver, &edata, *model);
    auto llh2 = rdata2->llh;
    ASSERT_NEAR(2.0 * llh1, llh2, 1e-6);

    // adjoint + replicates
    solver->set_sensitivity_method(amici::SensitivityMethod::adjoint);
    auto rdata3 = run_simulation(*solver, &edata, *model);
    auto llh3 = rdata3->llh;
    ASSERT_NEAR(llh2, llh3, 1e-6);
}
