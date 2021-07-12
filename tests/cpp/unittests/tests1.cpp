#include "testfunctions.h"

#include <amici/amici.h>
#include <amici/forwardproblem.h>
#include <amici/model_ode.h>
#include <amici/solver_cvodes.h>
#include <amici/solver_idas.h>
#include <amici/symbolic_functions.h>

#include <cmath>
#include <cstring>
#include <exception>
#include <vector>

#include <gtest/gtest.h>

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel()
{
    return std::make_unique<amici::Model_Test>();
}

} // namespace generic_model
} // namespace amici

using namespace amici;

namespace {

void
testSolverGetterSetters(CVodeSolver solver,
                        SensitivityMethod sensi_meth,
                        SensitivityOrder sensi,
                        InternalSensitivityMethod ism,
                        InterpolationType interp,
                        NonlinearSolverIteration iter,
                        LinearMultistepMethod lmm,
                        int steps,
                        int badsteps,
                        double tol,
                        double badtol);

class Model : public ::testing::Test {
  protected:
    // void SetUp() override {}
    // void TearDown() override {}

    int nx = 1, ny = 2, nz = 3, nmaxevent = 4;
    std::vector<realtype> p{ 1.0 };
    std::vector<realtype> k{ 0.5, 0.4, 0.7 };
    std::vector<int> plist{ 1 };
    std::vector<realtype> idlist{ 0 };
    std::vector<int> z2event{ 0, 0, 0 };
    Model_Test model = Model_Test(
        ModelDimensions(
            nx,        // nx_rdata
            nx,        // nxtrue_rdata
            nx,        // nx_solver
            nx,        // nxtrue_solver
            0,         // nx_solver_reinit
            static_cast<int>(p.size()),  // np
            static_cast<int>(k.size()),  // nk
            ny,        // ny
            ny,        // nytrue
            nz,        // nz
            nz,        // nztrue
            nmaxevent, // ne
            0,         // nJ
            0,         // nw
            0,         // ndwdx
            0,         // ndwdp
            0,         // dwdw
            0,         // ndxdotdw
            {},         // ndJydy
            0,         // nnz
            0,         // ubw
            0          // lbw
            ),
        SimulationParameters(k, p, plist),
        SecondOrderMode::none,
        idlist,
        z2event);
    std::vector<double> unscaled{ NAN };
};


TEST_F(Model, testScalingLin)
{
    model.setParameterScale(ParameterScaling::none);

    ASSERT_EQ(p[0], model.getParameters()[0]);
}

TEST_F(Model, testScalingLog)
{
    model.setParameterScale(ParameterScaling::ln);

    ASSERT_NEAR(std::log(p[0]), model.getParameters()[0], 1e-16);
}

TEST_F(Model, testScalingLog10)
{
    model.setParameterScale(ParameterScaling::log10);

    ASSERT_NEAR(std::log10(p[0]), model.getParameters()[0], 1e-16);
}

TEST_F(Model, testParameterScalingLengthMismatch)
{
    // too short
    auto pscale =
      std::vector<ParameterScaling>(p.size() - 1, ParameterScaling::log10);
    ASSERT_THROW(model.setParameterScale(pscale), AmiException);

    // too long
    pscale =
      std::vector<ParameterScaling>(p.size() + 1, ParameterScaling::log10);
    ASSERT_THROW(model.setParameterScale(pscale), AmiException);
}

TEST_F(Model, testSetTimepoints){
    ASSERT_THROW(model.setTimepoints(std::vector<realtype>{ 0.0, 1.0, 0.5 }),
                 AmiException);
}

TEST_F(Model, testNameIdGetterSetter)
{
    model.setParameterById("p0", 3.0);
    ASSERT_NEAR(model.getParameterById("p0"), 3.0, 1e-16);
    ASSERT_THROW(model.getParameterById("p1"), AmiException);
    ASSERT_NEAR(
      model.setParametersByIdRegex("p[\\d]+", 5.0), p.size(), 1e-16);
    for (const auto& ip : model.getParameters())
        ASSERT_NEAR(ip, 5.0, 1e-16);
    ASSERT_THROW(model.setParametersByIdRegex("k[\\d]+", 5.0), AmiException);

    model.setParameterByName("p0", 3.0);
    ASSERT_NEAR(model.getParameterByName("p0"), 3.0, 1e-16);
    ASSERT_THROW(model.getParameterByName("p1"), AmiException);
    ASSERT_NEAR(
      model.setParametersByNameRegex("p[\\d]+", 5.0), p.size(), 1e-16);
    for (const auto& ip : model.getParameters())
        ASSERT_NEAR(ip, 5.0, 1e-16);
    ASSERT_THROW(model.setParametersByNameRegex("k[\\d]+", 5.0), AmiException);

    model.setFixedParameterById("k0", 3.0);
    ASSERT_NEAR(model.getFixedParameterById("k0"), 3.0, 1e-16);
    ASSERT_THROW(model.getFixedParameterById("k4"), AmiException);
    ASSERT_NEAR(
      model.setFixedParametersByIdRegex("k[\\d]+", 5.0), k.size(), 1e-16);
    for (const auto& ik : model.getFixedParameters())
        ASSERT_NEAR(ik, 5.0, 1e-16);
    ASSERT_THROW(model.setFixedParametersByIdRegex("p[\\d]+", 5.0), AmiException);

    model.setFixedParameterByName("k0", 3.0);
    ASSERT_NEAR(model.getFixedParameterByName("k0"), 3.0, 1e-16);
    ASSERT_THROW(model.getFixedParameterByName("k4"), AmiException);
    ASSERT_NEAR(
      model.setFixedParametersByNameRegex("k[\\d]+", 5.0), k.size(), 1e-16);
    for (const auto& ik : model.getFixedParameters())
        ASSERT_NEAR(ik, 5.0, 1e-16);
    ASSERT_THROW(model.setFixedParametersByNameRegex("p[\\d]+", 5.0),
                 AmiException);
}

TEST_F(Model, reinitializeFixedParameterInitialStates)
{
    ASSERT_THROW(model.setReinitializeFixedParameterInitialStates(true),
                 AmiException);
    model.setReinitializeFixedParameterInitialStates(false);
    ASSERT_TRUE(!model.getReinitializeFixedParameterInitialStates());
    AmiVector x(nx);
    AmiVectorArray sx(model.np(), nx);
}

TEST(SymbolicFunctions, testSign)
{
    ASSERT_EQ(-1, sign(-2));
    ASSERT_EQ(0, sign(0));
    ASSERT_EQ(1, sign(2));
}

TEST(SymbolicFunctions, testHeaviside)
{
    ASSERT_EQ(0, heaviside(-1));
    ASSERT_EQ(1, heaviside(0));
    ASSERT_EQ(1, heaviside(1));
}

TEST(symbolicFunctions, testMin)
{
    ASSERT_EQ(-1, amici::min(-1, 2, 0));
    ASSERT_EQ(-2, amici::min(1, -2, 0));
    ASSERT_TRUE(amici::isNaN(amici::min(amici::getNaN(), amici::getNaN(), 0)));
    ASSERT_EQ(-1, amici::min(-1, amici::getNaN(), 0));
    ASSERT_EQ(-1, amici::min(amici::getNaN(), -1, 0));
}

TEST(symbolicFunctions, testMax)
{
    ASSERT_EQ(2, amici::max(-1, 2, 0));
    ASSERT_EQ(1, amici::max(1, -2, 0));
    ASSERT_TRUE(amici::isNaN(amici::max(amici::getNaN(), amici::getNaN(), 0)));
    ASSERT_EQ(-1, amici::max(-1, amici::getNaN(), 0));
    ASSERT_EQ(-1, amici::max(amici::getNaN(), -1, 0));
}

TEST(symbolicFunctions, testDMin)
{
    ASSERT_EQ(0, amici::Dmin(1, -1, -2, 0));
    ASSERT_EQ(1, amici::Dmin(1, -1, 2, 0));
    ASSERT_EQ(1, amici::Dmin(2, -1, -2, 0));
    ASSERT_EQ(0, amici::Dmin(2, -1, 2, 0));
}

TEST(symbolicFunctions, testDMax)
{
    ASSERT_EQ(1, amici::Dmax(1, -1, -2, 0));
    ASSERT_EQ(0, amici::Dmax(1, -1, 2, 0));
    ASSERT_EQ(0, amici::Dmax(2, -1, -2, 0));
    ASSERT_EQ(1, amici::Dmax(2, -1, 2, 0));
}

TEST(symbolicFunctions, testpos_pow)
{
    ASSERT_EQ(0, amici::pos_pow(-0.1, 3));
    ASSERT_EQ(pow(0.1, 3), amici::pos_pow(0.1, 3));
}

TEST(amiciSolver, testEquality)
{
    IDASolver i1, i2;
    CVodeSolver c1, c2;

    ASSERT_TRUE(i1 == i2);
    ASSERT_TRUE(c1 == c2);
    ASSERT_FALSE(i1 == c1);
}

TEST(amiciSolver, testClone)
{
    IDASolver i1;
    auto i2 = std::unique_ptr<Solver>(i1.clone());
    ASSERT_TRUE(i1 == *i2);

    CVodeSolver c1;
    auto c2 = std::unique_ptr<Solver>(c1.clone());
    ASSERT_TRUE(c1 == *c2);
    ASSERT_FALSE(*i2 == *c2);
}

TEST(amiciSolverIdas, testConstructionDestruction)
{
    IDASolver solver;
}

class edata : public ::testing::Test {
  protected:
    void SetUp() override {
        model->setTimepoints(timepoints);
        model->setNMaxEvent(nmaxevent);
        testModel.setTimepoints(timepoints);
        testModel.setNMaxEvent(nmaxevent);
    }

    int nx = 1, ny = 2, nz = 3, nmaxevent = 4;
    std::vector<realtype> timepoints = { 1, 2, 3, 4 };

    std::unique_ptr<amici::Model> model = amici::generic_model::getModel();

    Model_Test testModel = Model_Test(
                ModelDimensions(
                    nx,        // nx_rdata
                    nx,        // nxtrue_rdata
                    nx,        // nx_solver
                    nx,        // nxtrue_solver
                    0,         // nx_solver_reinit
                    1,  // np
                    3,  // nk
                    ny,        // ny
                    ny,        // nytrue
                    nz,        // nz
                    nz,        // nztrue
                    nmaxevent, // ne
                    0,         // nJ
                    0,         // nw
                    0,         // ndwdx
                    0,         // ndwdp
                    0,         // dwdw
                    0,         // ndxdotdw
                    {},         // ndJydy
                    0,         // nnz
                    0,         // ubw
                    0          // lbw
                    ),
                SimulationParameters(
                    std::vector<realtype>(3, 0.0),
                    std::vector<realtype>(1, 0.0),
                    std::vector<int>(2, 1)
                ),
                SecondOrderMode::none,
                std::vector<realtype>(),
                std::vector<int>());
};

TEST_F(edata, testConstructors1)
{
    auto edata = ExpData();
    ASSERT_TRUE(edata.nytrue() == 0);
    ASSERT_TRUE(edata.nztrue() == 0);
    ASSERT_TRUE(edata.nmaxevent() == 0);
}
TEST_F(edata, testConstructors2)
{
    auto edata = ExpData(model->nytrue, model->nztrue, model->nMaxEvent());
    ASSERT_TRUE(edata.nytrue() == model->nytrue);
    ASSERT_TRUE(edata.nztrue() == model->nztrue);
    ASSERT_TRUE(edata.nmaxevent() == model->nMaxEvent());
}

TEST_F(edata, testConstructors3)
{
    auto edata =
      ExpData(model->nytrue, model->nztrue, model->nMaxEvent(), timepoints);
    ASSERT_TRUE(edata.nytrue() == model->nytrue);
    ASSERT_TRUE(edata.nztrue() == model->nztrue);
    ASSERT_TRUE(edata.nmaxevent() == model->nMaxEvent());
    ASSERT_TRUE(edata.nt() == model->nt());
    checkEqualArray(
      timepoints, edata.getTimepoints(), TEST_ATOL, TEST_RTOL, "ts");
}

TEST_F(edata, testConstructors4)
{
    std::vector<realtype> y(ny * timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny * timepoints.size(), 0.1);
    std::vector<realtype> z(nz * nmaxevent, 0.0);
    std::vector<realtype> z_std(nz * nmaxevent, 0.1);

    auto edata = ExpData(testModel.nytrue,
                         testModel.nztrue,
                         testModel.nMaxEvent(),
                         timepoints,
                         y,
                         y_std,
                         z,
                         z_std);
    ASSERT_TRUE(edata.nytrue() == testModel.nytrue);
    ASSERT_TRUE(edata.nztrue() == testModel.nztrue);
    ASSERT_TRUE(edata.nmaxevent() == testModel.nMaxEvent());
    ASSERT_TRUE(edata.nt() == testModel.nt());
    checkEqualArray(
      timepoints, edata.getTimepoints(), TEST_ATOL, TEST_RTOL, "ts");
    checkEqualArray(
      y, edata.getObservedData(), TEST_ATOL, TEST_RTOL, "observedData");
    checkEqualArray(y_std,
                    edata.getObservedDataStdDev(),
                    TEST_ATOL,
                    TEST_RTOL,
                    "observedDataStdDev");
    checkEqualArray(
      z, edata.getObservedEvents(), TEST_ATOL, TEST_RTOL, "observedEvents");
    checkEqualArray(z_std,
                    edata.getObservedEventsStdDev(),
                    TEST_ATOL,
                    TEST_RTOL,
                    "observedEventsStdDev");

    auto edata_copy = ExpData(edata);
    ASSERT_TRUE(edata.nytrue() == edata_copy.nytrue());
    ASSERT_TRUE(edata.nztrue() == edata_copy.nztrue());
    ASSERT_TRUE(edata.nmaxevent() == edata_copy.nmaxevent());
    ASSERT_TRUE(edata.nt() == edata_copy.nt());
    checkEqualArray(edata_copy.getTimepoints(),
                    edata.getTimepoints(),
                    TEST_ATOL,
                    TEST_RTOL,
                    "ts");
    checkEqualArray(edata_copy.getObservedData(),
                    edata.getObservedData(),
                    TEST_ATOL,
                    TEST_RTOL,
                    "observedData");
    checkEqualArray(edata_copy.getObservedDataStdDev(),
                    edata.getObservedDataStdDev(),
                    TEST_ATOL,
                    TEST_RTOL,
                    "observedDataStdDev");
    checkEqualArray(edata_copy.getObservedEvents(),
                    edata.getObservedEvents(),
                    TEST_ATOL,
                    TEST_RTOL,
                    "observedEvents");
    checkEqualArray(edata_copy.getObservedEventsStdDev(),
                    edata.getObservedEventsStdDev(),
                    TEST_ATOL,
                    TEST_RTOL,
                    "observedEventsStdDev");
}

TEST_F(edata, testConstructors5)
{
    testModel.setTimepoints(timepoints);
    auto edata = ExpData(testModel);
    ASSERT_TRUE(edata.nytrue() == testModel.nytrue);
    ASSERT_TRUE(edata.nztrue() == testModel.nztrue);
    ASSERT_TRUE(edata.nmaxevent() == testModel.nMaxEvent());
    ASSERT_TRUE(edata.nt() == testModel.nt());
    checkEqualArray(testModel.getTimepoints(),
                    edata.getTimepoints(),
                    TEST_ATOL,
                    TEST_RTOL,
                    "ts");
}

TEST_F(edata, testDimensionChecks)
{

    std::vector<realtype> bad_std(ny, -0.1);

    std::vector<realtype> y(ny * timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny * timepoints.size(), 0.1);
    std::vector<realtype> z(nz * nmaxevent, 0.0);
    std::vector<realtype> z_std(nz * nmaxevent, 0.1);

    ASSERT_THROW(ExpData(testModel.nytrue,
                         testModel.nztrue,
                         testModel.nMaxEvent(),
                         timepoints,
                         z,
                         z_std,
                         z,
                         z_std),
                 AmiException);

    ASSERT_THROW(ExpData(testModel.nytrue,
                         testModel.nztrue,
                         testModel.nMaxEvent(),
                         timepoints,
                         z,
                         bad_std,
                         z,
                         z_std),
                 AmiException);

    auto edata = ExpData(testModel);

    std::vector<realtype> bad_y(ny * timepoints.size() + 1, 0.0);
    std::vector<realtype> bad_y_std(ny * timepoints.size() + 1, 0.1);
    std::vector<realtype> bad_z(nz * nmaxevent + 1, 0.0);
    std::vector<realtype> bad_z_std(nz * nmaxevent + 1, 0.1);

    ASSERT_THROW(edata.setObservedData(bad_y), AmiException);
    ASSERT_THROW(edata.setObservedDataStdDev(bad_y_std), AmiException);
    ASSERT_THROW(edata.setObservedEvents(bad_z), AmiException);
    ASSERT_THROW(edata.setObservedEventsStdDev(bad_y_std), AmiException);

    std::vector<realtype> bad_single_y(edata.nt() + 1, 0.0);
    std::vector<realtype> bad_single_y_std(edata.nt() + 1, 0.1);
    std::vector<realtype> bad_single_z(edata.nmaxevent() + 1, 0.0);
    std::vector<realtype> bad_single_z_std(edata.nmaxevent() + 1, 0.1);

    ASSERT_THROW(edata.setObservedData(bad_single_y, 0),
                 AmiException);
    ASSERT_THROW(edata.setObservedDataStdDev(bad_single_y_std, 0),
                 AmiException);
    ASSERT_THROW(edata.setObservedEvents(bad_single_z, 0),
                 AmiException);
    ASSERT_THROW(edata.setObservedEventsStdDev(bad_single_y_std, 0),
                 AmiException);

    ASSERT_THROW(edata.setTimepoints(std::vector<realtype>{ 0.0, 1.0, 0.5 }),
                 AmiException);
}

TEST_F(edata, testSettersGetters)
{
    auto edata = ExpData(testModel);

    std::vector<realtype> y(ny * timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny * timepoints.size(), 0.1);
    std::vector<realtype> z(nz * nmaxevent, 0.0);
    std::vector<realtype> z_std(nz * nmaxevent, 0.1);

    edata.setObservedData(y);
    checkEqualArray(
      edata.getObservedData(), y, TEST_ATOL, TEST_RTOL, "ObservedData");
    edata.setObservedDataStdDev(y_std);
    checkEqualArray(edata.getObservedDataStdDev(),
                    y_std,
                    TEST_ATOL,
                    TEST_RTOL,
                    "ObservedDataStdDev");
    edata.setObservedEvents(z);
    checkEqualArray(
      edata.getObservedEvents(), z, TEST_ATOL, TEST_RTOL, "ObservedEvents");
    edata.setObservedEventsStdDev(z_std);
    checkEqualArray(edata.getObservedEventsStdDev(),
                    z_std,
                    TEST_ATOL,
                    TEST_RTOL,
                    "ObservedEventsStdDev");

    std::vector<realtype> single_y(edata.nt(), 0.0);
    std::vector<realtype> single_y_std(edata.nt(), 0.1);

    for (int iy = 0; iy < ny; ++iy) {
        edata.setObservedData(single_y, iy);
        edata.setObservedDataStdDev(single_y_std, iy);
    }
    ASSERT_THROW(edata.setObservedData(single_y, ny), std::exception);
    ASSERT_THROW(edata.setObservedData(single_y, -1), std::exception);
    ASSERT_THROW(edata.setObservedDataStdDev(single_y_std, ny), std::exception);
    ASSERT_THROW(edata.setObservedDataStdDev(single_y_std, -1), std::exception);

    std::vector<realtype> single_z(edata.nmaxevent(), 0.0);
    std::vector<realtype> single_z_std(edata.nmaxevent(), 0.1);

    for (int iz = 0; iz < nz; ++iz) {
        edata.setObservedEvents(single_z, iz);
        edata.setObservedEventsStdDev(single_z_std, iz);
    }

    ASSERT_THROW(edata.setObservedEvents(single_z, nz), std::exception);
    ASSERT_THROW(edata.setObservedEvents(single_z, -1), std::exception);
    ASSERT_THROW(edata.setObservedEventsStdDev(single_z_std, nz),
                 std::exception);
    ASSERT_THROW(edata.setObservedEventsStdDev(single_z_std, -1),
                 std::exception);

    ASSERT_TRUE(edata.getObservedDataPtr(0));
    ASSERT_TRUE(edata.getObservedDataStdDevPtr(0));
    ASSERT_TRUE(edata.getObservedEventsPtr(0));
    ASSERT_TRUE(edata.getObservedEventsStdDevPtr(0));

    std::vector<realtype> empty(0, 0.0);

    edata.setObservedData(empty);
    edata.setObservedDataStdDev(empty);
    edata.setObservedEvents(empty);
    edata.setObservedEventsStdDev(empty);

    ASSERT_TRUE(!edata.getObservedDataPtr(0));
    ASSERT_TRUE(!edata.getObservedDataStdDevPtr(0));
    ASSERT_TRUE(!edata.getObservedEventsPtr(0));
    ASSERT_TRUE(!edata.getObservedEventsStdDevPtr(0));

    checkEqualArray(
      edata.getObservedData(), empty, TEST_ATOL, TEST_RTOL, "ObservedData");
    checkEqualArray(edata.getObservedDataStdDev(),
                    empty,
                    TEST_ATOL,
                    TEST_RTOL,
                    "ObservedDataStdDev");
    checkEqualArray(
      edata.getObservedEvents(), empty, TEST_ATOL, TEST_RTOL, "ObservedEvents");
    checkEqualArray(edata.getObservedEventsStdDev(),
                    empty,
                    TEST_ATOL,
                    TEST_RTOL,
                    "ObservedEventsStdDev");
}

class SolverTest : public ::testing::Test {
  protected:
    void SetUp() override {
        tol = 0.01;
        badtol = -0.01;
        sensi_meth = SensitivityMethod::adjoint;
        sensi = SensitivityOrder::first;
        steps = 1000;
        badsteps = -1;
        lmm = LinearMultistepMethod::adams;
        iter = NonlinearSolverIteration::fixedpoint;
        ism = InternalSensitivityMethod::staggered1;
        interp = InterpolationType::polynomial;
    }

    int nx = 1, ny = 2, nz = 3, ne = 0;
    double tol, badtol;
    std::vector<realtype> timepoints = { 1, 2, 3, 4 };

    std::unique_ptr<amici::Model> model = amici::generic_model::getModel();
    SensitivityMethod sensi_meth;
    SensitivityOrder sensi;
    int steps, badsteps;
    LinearMultistepMethod lmm;
    NonlinearSolverIteration iter;
    InternalSensitivityMethod ism;
    InterpolationType interp;

    Model_Test testModel = Model_Test(
                ModelDimensions(
                    nx,        // nx_rdata
                    nx,        // nxtrue_rdata
                    nx,        // nx_solver
                    nx,        // nxtrue_solver
                    0,         // nx_solver_reinit
                    1,         // np
                    3,         // nk
                    ny,        // ny
                    ny,        // nytrue
                    nz,        // nz
                    nz,        // nztrue
                    ne,        // ne
                    0,         // nJ
                    0,         // nw
                    0,         // ndwdx
                    0,         // ndwdp
                    0,         // dwdw
                    0,         // ndxdotdw
                    {},         // ndJydy
                    1,         // nnz
                    0,         // ubw
                    0         // lbw
                    ),
                SimulationParameters(
                    std::vector<realtype>(3, 0.0),
                    std::vector<realtype>(1, 0.0),
                    std::vector<int>(2, 1)
                ),
                SecondOrderMode::none,
                std::vector<realtype>(0, 0.0),
                std::vector<int>());

            CVodeSolver solver = CVodeSolver();
};

TEST_F(SolverTest, testSettersGettersNoSetup)
{
    testSolverGetterSetters(solver,
                            sensi_meth,
                            sensi,
                            ism,
                            interp,
                            iter,
                            lmm,
                            steps,
                            badsteps,
                            tol,
                            badtol);
}

TEST_F(SolverTest, testSettersGettersWithSetup)
{

    solver.setSensitivityMethod(sensi_meth);
    ASSERT_EQ(static_cast<int>(solver.getSensitivityMethod()),
                static_cast<int>(sensi_meth));

    auto rdata =
      std::unique_ptr<ReturnData>(new ReturnData(solver, testModel));
    AmiVector x(nx), dx(nx);
    AmiVectorArray sx(nx, 1), sdx(nx, 1);

    testModel.setInitialStates(std::vector<realtype>{ 0 });

    solver.setup(0, &testModel, x, dx, sx, sdx);

    testSolverGetterSetters(solver,
                            sensi_meth,
                            sensi,
                            ism,
                            interp,
                            iter,
                            lmm,
                            steps,
                            badsteps,
                            tol,
                            badtol);
}

void
testSolverGetterSetters(CVodeSolver solver,
                        SensitivityMethod sensi_meth,
                        SensitivityOrder sensi,
                        InternalSensitivityMethod ism,
                        InterpolationType interp,
                        NonlinearSolverIteration iter,
                        LinearMultistepMethod lmm,
                        int steps,
                        int badsteps,
                        double tol,
                        double badtol)
{

    solver.setSensitivityMethod(sensi_meth);
    ASSERT_EQ(static_cast<int>(solver.getSensitivityMethod()),
                static_cast<int>(sensi_meth));

    solver.setSensitivityOrder(sensi);
    ASSERT_EQ(static_cast<int>(solver.getSensitivityOrder()),
                static_cast<int>(sensi));

    solver.setInternalSensitivityMethod(ism);
    ASSERT_EQ(static_cast<int>(solver.getInternalSensitivityMethod()),
                static_cast<int>(ism));

    solver.setInterpolationType(interp);
    ASSERT_EQ(static_cast<int>(solver.getInterpolationType()),
                static_cast<int>(interp));

    solver.setNonlinearSolverIteration(iter);
    ASSERT_EQ(static_cast<int>(solver.getNonlinearSolverIteration()),
                static_cast<int>(iter));

    solver.setLinearMultistepMethod(lmm);
    ASSERT_EQ(static_cast<int>(solver.getLinearMultistepMethod()),
                static_cast<int>(lmm));

    solver.setPreequilibration(true);
    ASSERT_EQ(solver.getPreequilibration(), true);

    solver.setStabilityLimitFlag(true);
    ASSERT_EQ(solver.getStabilityLimitFlag(), true);

    ASSERT_THROW(solver.setNewtonMaxSteps(badsteps), AmiException);
    solver.setNewtonMaxSteps(steps);
    ASSERT_EQ(solver.getNewtonMaxSteps(), steps);

    ASSERT_THROW(solver.setNewtonMaxLinearSteps(badsteps), AmiException);
    solver.setNewtonMaxLinearSteps(steps);
    ASSERT_EQ(solver.getNewtonMaxLinearSteps(), steps);

    ASSERT_THROW(solver.setMaxSteps(badsteps), AmiException);
    solver.setMaxSteps(steps);
    ASSERT_EQ(solver.getMaxSteps(), steps);

    ASSERT_THROW(solver.setMaxStepsBackwardProblem(badsteps), AmiException);
    solver.setMaxStepsBackwardProblem(steps);
    ASSERT_EQ(solver.getMaxStepsBackwardProblem(), steps);

    ASSERT_THROW(solver.setRelativeTolerance(badtol), AmiException);
    solver.setRelativeTolerance(tol);
    ASSERT_EQ(solver.getRelativeTolerance(), tol);

    ASSERT_THROW(solver.setAbsoluteTolerance(badtol), AmiException);
    solver.setAbsoluteTolerance(tol);
    ASSERT_EQ(solver.getAbsoluteTolerance(), tol);

    ASSERT_THROW(solver.setRelativeToleranceQuadratures(badtol), AmiException);
    solver.setRelativeToleranceQuadratures(tol);
    ASSERT_EQ(solver.getRelativeToleranceQuadratures(), tol);

    ASSERT_THROW(solver.setAbsoluteToleranceQuadratures(badtol), AmiException);
    solver.setAbsoluteToleranceQuadratures(tol);
    ASSERT_EQ(solver.getAbsoluteToleranceQuadratures(), tol);

    ASSERT_THROW(solver.setRelativeToleranceSteadyState(badtol), AmiException);
    solver.setRelativeToleranceSteadyState(tol);
    ASSERT_EQ(solver.getRelativeToleranceSteadyState(), tol);

    ASSERT_THROW(solver.setAbsoluteToleranceSteadyState(badtol), AmiException);
    solver.setAbsoluteToleranceSteadyState(tol);
    ASSERT_EQ(solver.getAbsoluteToleranceSteadyState(), tol);
}

class AmiVectorTest : public ::testing::Test {
  protected:
    std::vector<double> vec1{ 1, 2, 4, 3 };
    std::vector<double> vec2{ 4, 1, 2, 3 };
    std::vector<double> vec3{ 4, 4, 2, 1 };
};

TEST_F(AmiVectorTest, vector)
{
    AmiVector av(vec1);
    N_Vector nvec = av.getNVector();
    for (int i = 0; i < av.getLength(); ++i)
        ASSERT_EQ(av.at(i), NV_Ith_S(nvec, i));
}

TEST_F(AmiVectorTest, vectorArray)
{
    AmiVectorArray ava(4, 3);
    AmiVector av1(vec1), av2(vec2), av3(vec3);
    std::vector<AmiVector> avs{ av1, av2, av3 };
    for (int i = 0; i < ava.getLength(); ++i)
        ava[i] = avs.at(i);

    std::vector<double> badLengthVector(13, 0.0);
    std::vector<double> flattened(12, 0.0);

    ASSERT_THROW(ava.flatten_to_vector(badLengthVector), AmiException);
    ava.flatten_to_vector(flattened);
    for (int i = 0; i < ava.getLength(); ++i) {
        const AmiVector av = ava[i];
        for (int j = 0; j < av.getLength(); ++j)
            ASSERT_EQ(flattened.at(i * av.getLength() + j), av.at(j));
    }
}

class SunMatrixWrapperTest : public ::testing::Test {
  protected:
    void SetUp() override {
        A.set_data(0, 0, 0.69);
        A.set_data(1, 0, 0.32);
        A.set_data(2, 0, 0.95);
        A.set_data(0, 1, 0.03);
        A.set_data(1, 1, 0.44);
        A.set_data(2, 1, 0.38);

        B.set_indexptr(0, 0);
        B.set_indexptr(1, 2);
        B.set_indexptr(2, 4);
        B.set_indexptr(3, 5);
        B.set_indexptr(4, 7);
        B.set_data(0, 3);
        B.set_data(1, 1);
        B.set_data(2, 3);
        B.set_data(3, 7);
        B.set_data(4, 1);
        B.set_data(5, 2);
        B.set_data(6, 9);
        B.set_indexval(0, 1);
        B.set_indexval(1, 3);
        B.set_indexval(2, 0);
        B.set_indexval(3, 2);
        B.set_indexval(4, 0);
        B.set_indexval(5, 1);
        B.set_indexval(6, 3);
    }

    //inputs
    std::vector<double> a{0.82, 0.91, 0.13};
    std::vector<double> b{0.77, 0.80};
    SUNMatrixWrapper A = SUNMatrixWrapper(3, 2);
    SUNMatrixWrapper B = SUNMatrixWrapper(4, 4, 7, CSC_MAT);
    // result
    std::vector<double> d{1.3753, 1.5084, 1.1655};
};

TEST_F(SunMatrixWrapperTest, sparse_multiply)
{

    auto A_sparse = SUNMatrixWrapper(A, 0.0, CSC_MAT);
    auto c(a); //copy c
    A_sparse.multiply(c, b);
    checkEqualArray(d, c, TEST_ATOL, TEST_RTOL, "multiply");
}

TEST_F(SunMatrixWrapperTest, sparse_multiply_empty)
{
    // Ensure empty Matrix vector multiplication succeeds
    auto A_sparse = SUNMatrixWrapper(1, 1, 0, CSR_MAT);
    std::vector<double> b {0.1};
    std::vector<double> c {0.1};
    A_sparse.multiply(c, b);
    ASSERT_TRUE(c[0] == 0.1);

    A_sparse = SUNMatrixWrapper(1, 1, 0, CSC_MAT);
    A_sparse.multiply(c, b);
    ASSERT_TRUE(c[0] == 0.1);
}

TEST_F(SunMatrixWrapperTest, dense_multiply)
{
    auto c(a); //copy c
    A.multiply(c, b);
    checkEqualArray(d, c, TEST_ATOL, TEST_RTOL, "multiply");
}

TEST_F(SunMatrixWrapperTest, multiply_throws)
{
    auto b_amivector = AmiVector(b);
    auto a_amivector = AmiVector(a);
}

TEST_F(SunMatrixWrapperTest, transform_throws)
{
    ASSERT_THROW(SUNMatrixWrapper(A, 0.0, 13), std::invalid_argument);
    auto A_sparse = SUNMatrixWrapper(A, 0.0, CSR_MAT);
    ASSERT_THROW(SUNMatrixWrapper(A_sparse, 0.0, CSR_MAT),
                 std::invalid_argument);
}

TEST_F(SunMatrixWrapperTest, block_transpose)
{
    auto B_sparse = SUNMatrixWrapper(4, 4, 7, CSR_MAT);
    ASSERT_THROW(B.transpose(B_sparse, 1.0, 4), std::domain_error);

    B_sparse = SUNMatrixWrapper(4, 4, 7, CSC_MAT);
    B.transpose(B_sparse, -1.0, 2);
    for (int idx = 0; idx < 7; idx++) {
        ASSERT_TRUE(SM_INDEXVALS_S(B.get())[idx]
                   == SM_INDEXVALS_S(B_sparse.get())[idx]);
        if (idx == 1) {
            ASSERT_TRUE(SM_DATA_S(B.get())[idx]
                       == -SM_DATA_S(B_sparse.get())[3]);
        } else if (idx == 3) {
            ASSERT_TRUE(SM_DATA_S(B.get())[idx]
                       == -SM_DATA_S(B_sparse.get())[1]);
        } else {
            ASSERT_TRUE(SM_DATA_S(B.get())[idx]
                       == -SM_DATA_S(B_sparse.get())[idx]);
        }
    }
    for (int icol = 0; icol <= 4; icol++)
        ASSERT_TRUE(SM_INDEXPTRS_S(B.get())[icol]
                   == SM_INDEXPTRS_S(B_sparse.get())[icol]);
}

} // namespace
