#include "testfunctions.h"

#include <amici/amici.h>
#include <amici/model_ode.h>
#include <amici/solver_cvodes.h>
#include <amici/symbolic_functions.h>

#include <exception>
#include <vector>

#include <gtest/gtest.h>

namespace amici {
namespace generic_model {
std::unique_ptr<Model> get_model();
} // namespace generic_model
} // namespace amici

using namespace amici;

namespace {

class ExpDataTest : public ::testing::Test {
  protected:
    void SetUp() override {
        model->set_timepoints(timepoints);
        model->set_n_max_event(nmaxevent);
        testModel.set_timepoints(timepoints);
        testModel.set_n_max_event(nmaxevent);
    }

    int nx = 1, ny = 2, nz = 3, nmaxevent = 4;
    std::vector<realtype> timepoints = {1, 2, 3, 4};

    std::unique_ptr<Model> model = generic_model::get_model();

    Model_Test testModel = Model_Test(
        ModelDimensions{
            .nx_rdata = nx,
            .nxtrue_rdata = nx,
            .nx_solver = nx,
            .nxtrue_solver = nx,
            .nx_solver_reinit = 0,
            .np = 1,
            .nk = 3,
            .ny = ny,
            .nytrue = ny,
            .nz = nz,
            .nztrue = nz,
            .ne = nmaxevent,
            .ne_solver = 0,
            .nspl = 0,
            .nw = 0,
            .ndwdx = 0,
            .ndwdp = 0,
            .ndwdw = 0,
            .ndxdotdw = 0,
            .ndJydy = {0, 0},
            .ndxrdatadxsolver = 0,
            .ndxrdatadtcl = 0,
            .ndtotal_cldx_rdata = 0,
            .nnz = 0,
            .nJ = 0,
            .ubw = 0,
            .lbw = 0
        },
        SimulationParameters(
            std::vector<realtype>(3, 0.0), std::vector<realtype>(1, 0.0),
            std::vector<int>(2, 1)
        ),
        SecondOrderMode::none, std::vector<realtype>(), std::vector<int>(),
        {
            Event("e1", true, true, 1),
            Event("e1", true, true, 1),
            Event("e1", true, true, 1),
            Event("e1", true, true, 1),
        }
    );
};

TEST_F(ExpDataTest, DefaultConstructable) {
    ExpData edata{};
    ASSERT_EQ(edata.nytrue(), 0);
    ASSERT_EQ(edata.nztrue(), 0);
    ASSERT_EQ(edata.nmaxevent(), 0);
}
TEST_F(ExpDataTest, ModelCtor) {
    ExpData edata(model->nytrue, model->nztrue, model->n_max_event());
    ASSERT_EQ(edata.nytrue(), model->nytrue);
    ASSERT_EQ(edata.nztrue(), model->nztrue);
    ASSERT_EQ(edata.nmaxevent(), model->n_max_event());
}

TEST_F(ExpDataTest, DimensionCtor) {
    ExpData edata(
        model->nytrue, model->nztrue, model->n_max_event(), timepoints
    );
    ASSERT_EQ(edata.nytrue(), model->nytrue);
    ASSERT_EQ(edata.nztrue(), model->nztrue);
    ASSERT_EQ(edata.nmaxevent(), model->n_max_event());
    ASSERT_EQ(edata.nt(), model->nt());
    checkEqualArray(
        timepoints, edata.get_timepoints(), TEST_ATOL, TEST_RTOL, "ts"
    );
}

TEST_F(ExpDataTest, MeasurementCtor) {
    std::vector<realtype> y(ny * timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny * timepoints.size(), 0.1);
    std::vector<realtype> z(nz * nmaxevent, 0.0);
    std::vector<realtype> z_std(nz * nmaxevent, 0.1);

    ExpData edata(
        testModel.nytrue, testModel.nztrue, testModel.n_max_event(), timepoints,
        y, y_std, z, z_std
    );
    ASSERT_EQ(edata.nytrue(), testModel.nytrue);
    ASSERT_EQ(edata.nztrue(), testModel.nztrue);
    ASSERT_EQ(edata.nmaxevent(), testModel.n_max_event());
    ASSERT_EQ(edata.nt(), testModel.nt());
    checkEqualArray(
        timepoints, edata.get_timepoints(), TEST_ATOL, TEST_RTOL, "ts"
    );
    checkEqualArray(
        y, edata.get_measurements(), TEST_ATOL, TEST_RTOL, "observedData"
    );
    checkEqualArray(
        y_std, edata.get_noise_scales(), TEST_ATOL, TEST_RTOL,
        "observedDataStdDev"
    );
    checkEqualArray(
        z, edata.get_event_measurements(), TEST_ATOL, TEST_RTOL, "observedEvents"
    );
    checkEqualArray(
        z_std, edata.get_event_noise_scales(), TEST_ATOL, TEST_RTOL,
        "observedEventsStdDev"
    );

    ExpData edata_copy(edata);
    ASSERT_EQ(edata.nytrue(), edata_copy.nytrue());
    ASSERT_EQ(edata.nztrue(), edata_copy.nztrue());
    ASSERT_EQ(edata.nmaxevent(), edata_copy.nmaxevent());
    ASSERT_EQ(edata.nt(), edata_copy.nt());
    checkEqualArray(
        edata_copy.get_timepoints(), edata.get_timepoints(), TEST_ATOL,
        TEST_RTOL, "ts"
    );
    checkEqualArray(
        edata_copy.get_measurements(), edata.get_measurements(), TEST_ATOL,
        TEST_RTOL, "observedData"
    );
    checkEqualArray(
        edata_copy.get_noise_scales(),
        edata.get_noise_scales(), TEST_ATOL, TEST_RTOL,
        "observedDataStdDev"
    );
    checkEqualArray(
        edata_copy.get_event_measurements(), edata.get_event_measurements(),
        TEST_ATOL, TEST_RTOL, "observedEvents"
    );
    checkEqualArray(
        edata_copy.get_event_noise_scales(),
        edata.get_event_noise_scales(), TEST_ATOL, TEST_RTOL,
        "observedEventsStdDev"
    );
}

TEST_F(ExpDataTest, CopyConstructable) {
    testModel.set_timepoints(timepoints);
    auto edata = ExpData(testModel);
    ASSERT_EQ(edata.nytrue(), testModel.nytrue);
    ASSERT_EQ(edata.nztrue(), testModel.nztrue);
    ASSERT_EQ(edata.nmaxevent(), testModel.n_max_event());
    ASSERT_EQ(edata.nt(), testModel.nt());
    checkEqualArray(
        testModel.get_timepoints(), edata.get_timepoints(), TEST_ATOL,
        TEST_RTOL, "ts"
    );
}

TEST_F(ExpDataTest, Equality) {
    auto edata = ExpData(testModel);
    auto edata2(edata);
    ASSERT_TRUE(edata == edata2);

    edata2.id = "different";
    ASSERT_FALSE(edata == edata2);
}

TEST_F(ExpDataTest, DimensionChecks) {
    std::vector<realtype> bad_std(ny, -0.1);
    std::vector<realtype> y(ny * timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny * timepoints.size(), 0.1);
    std::vector<realtype> z(nz * nmaxevent, 0.0);
    std::vector<realtype> z_std(nz * nmaxevent, 0.1);

    ASSERT_THROW(
        ExpData(
            testModel.nytrue, testModel.nztrue, testModel.n_max_event(),
            timepoints, z, z_std, z, z_std
        ),
        AmiException
    );

    ASSERT_THROW(
        ExpData(
            testModel.nytrue, testModel.nztrue, testModel.n_max_event(),
            timepoints, z, bad_std, z, z_std
        ),
        AmiException
    );

    ExpData edata(testModel);

    std::vector<realtype> bad_y(ny * timepoints.size() + 1, 0.0);
    std::vector<realtype> bad_y_std(ny * timepoints.size() + 1, 0.1);
    std::vector<realtype> bad_z(nz * nmaxevent + 1, 0.0);
    std::vector<realtype> bad_z_std(nz * nmaxevent + 1, 0.1);

    ASSERT_THROW(edata.set_measurements(bad_y), AmiException);
    ASSERT_THROW(edata.set_noise_scales(bad_y_std), AmiException);
    ASSERT_THROW(edata.set_event_measurements(bad_z), AmiException);
    ASSERT_THROW(edata.set_event_noise_scales(bad_y_std), AmiException);

    std::vector<realtype> bad_single_y(edata.nt() + 1, 0.0);
    std::vector<realtype> bad_single_y_std(edata.nt() + 1, 0.1);
    std::vector<realtype> bad_single_z(edata.nmaxevent() + 1, 0.0);
    std::vector<realtype> bad_single_z_std(edata.nmaxevent() + 1, 0.1);

    ASSERT_THROW(edata.set_measurements(bad_single_y, 0), AmiException);
    ASSERT_THROW(
        edata.set_noise_scales(bad_single_y_std, 0), AmiException
    );
    ASSERT_THROW(edata.set_event_measurements(bad_single_z, 0), AmiException);
    ASSERT_THROW(
        edata.set_event_noise_scales(bad_single_y_std, 0), AmiException
    );

    ASSERT_THROW(
        edata.set_timepoints(std::vector<realtype>{0.0, 1.0, 0.5}), AmiException
    );
}

TEST_F(ExpDataTest, SettersGetters) {
    ExpData edata(testModel);

    std::vector<realtype> y(ny * timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny * timepoints.size(), 0.1);
    std::vector<realtype> z(nz * nmaxevent, 0.0);
    std::vector<realtype> z_std(nz * nmaxevent, 0.1);

    edata.set_measurements(y);
    checkEqualArray(
        edata.get_measurements(), y, TEST_ATOL, TEST_RTOL, "ObservedData"
    );
    edata.set_noise_scales(y_std);
    checkEqualArray(
        edata.get_noise_scales(), y_std, TEST_ATOL, TEST_RTOL,
        "ObservedDataStdDev"
    );
    edata.set_event_measurements(z);
    checkEqualArray(
        edata.get_event_measurements(), z, TEST_ATOL, TEST_RTOL, "ObservedEvents"
    );
    edata.set_event_noise_scales(z_std);
    checkEqualArray(
        edata.get_event_noise_scales(), z_std, TEST_ATOL, TEST_RTOL,
        "ObservedEventsStdDev"
    );

    std::vector<realtype> single_y(edata.nt(), 0.0);
    std::vector<realtype> single_y_std(edata.nt(), 0.1);

    for (int iy = 0; iy < ny; ++iy) {
        edata.set_measurements(single_y, iy);
        edata.set_noise_scales(single_y_std, iy);
    }
    ASSERT_THROW(edata.set_measurements(single_y, ny), std::exception);
    ASSERT_THROW(edata.set_measurements(single_y, -1), std::exception);
    ASSERT_THROW(
        edata.set_noise_scales(single_y_std, ny), std::exception
    );
    ASSERT_THROW(
        edata.set_noise_scales(single_y_std, -1), std::exception
    );

    std::vector<realtype> single_z(edata.nmaxevent(), 0.0);
    std::vector<realtype> single_z_std(edata.nmaxevent(), 0.1);

    for (int iz = 0; iz < nz; ++iz) {
        edata.set_event_measurements(single_z, iz);
        edata.set_event_noise_scales(single_z_std, iz);
    }

    ASSERT_THROW(edata.set_event_measurements(single_z, nz), std::exception);
    ASSERT_THROW(edata.set_event_measurements(single_z, -1), std::exception);
    ASSERT_THROW(
        edata.set_event_noise_scales(single_z_std, nz), std::exception
    );
    ASSERT_THROW(
        edata.set_event_noise_scales(single_z_std, -1), std::exception
    );

    ASSERT_TRUE(edata.get_measurements_ptr(0));
    ASSERT_TRUE(edata.get_noise_scales_ptr(0));
    ASSERT_TRUE(edata.get_event_measurements_ptr(0));
    ASSERT_TRUE(edata.get_event_noise_scales_ptr(0));

    std::vector<realtype> empty(0, 0.0);

    edata.set_measurements(empty);
    edata.set_noise_scales(empty);
    edata.set_event_measurements(empty);
    edata.set_event_noise_scales(empty);

    ASSERT_TRUE(!edata.get_measurements_ptr(0));
    ASSERT_TRUE(!edata.get_noise_scales_ptr(0));
    ASSERT_TRUE(!edata.get_event_measurements_ptr(0));
    ASSERT_TRUE(!edata.get_event_noise_scales_ptr(0));

    checkEqualArray(
        edata.get_measurements(), empty, TEST_ATOL, TEST_RTOL, "ObservedData"
    );
    checkEqualArray(
        edata.get_noise_scales(), empty, TEST_ATOL, TEST_RTOL,
        "ObservedDataStdDev"
    );
    checkEqualArray(
        edata.get_event_measurements(), empty, TEST_ATOL, TEST_RTOL,
        "ObservedEvents"
    );
    checkEqualArray(
        edata.get_event_noise_scales(), empty, TEST_ATOL, TEST_RTOL,
        "ObservedEventsStdDev"
    );
}

TEST_F(ExpDataTest, RngSeed) {
    ReturnData rdata(CVodeSolver(), testModel);
    rdata.y.assign(testModel.ny * testModel.nt(), 1.0);
    rdata.z.assign(testModel.nz * testModel.n_max_event(), 1.0);

    // random RNG seed
    ExpData edata1(rdata, 1, 1, -1);

    ASSERT_TRUE(edata1.get_measurements() != rdata.y);

    // fixed RNG seed
    ExpData edata2(rdata, 1, 1, 1);

    ASSERT_TRUE(edata2.get_measurements() != rdata.y);
    ASSERT_TRUE(edata2.get_measurements() != edata1.get_measurements());

    ExpData edata3(rdata, 1, 1, 1);
    ASSERT_TRUE(edata3.get_measurements() == edata2.get_measurements());
}

} // namespace
