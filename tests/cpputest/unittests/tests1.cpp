#include "testfunctions.h"

#include <amici/amici.h>
#include <amici/solver_idas.h>
#include <amici/solver_cvodes.h>
#include <amici/symbolic_functions.h>
#include <amici/model_ode.h>
#include <amici/forwardproblem.h>

#include <cstring>
#include <cmath>
#include <vector>
#include <exception>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(new amici::Model_Test());
}

using namespace amici;

void testSolverGetterSetters(CVodeSolver solver,
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

TEST_GROUP(amici)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST_GROUP(model)
{
    int nx = 1, ny = 2, nz = 3, nmaxevent = 4;
    std::vector<realtype> p {1.0};
    std::vector<realtype> k {0.5,0.4,0.7};
    std::vector<int> plist {1};
    std::vector<realtype> idlist {0};
    std::vector<int> z2event {0,0,0};
    Model_Test model = Model_Test(nx, nx, nx, nx, ny, ny, nz, nz, nmaxevent,
                                  0, 0, 0, 0, 0, 0, 0, SecondOrderMode::none,
                                  p, k, plist, idlist, z2event);
    
    std::vector<double> unscaled {NAN};
    
    void setup() {
    }

    void teardown() {

    }
};

TEST(model, testScalingLin) {
    model.setParameterScale(ParameterScaling::none);

    CHECK_EQUAL(p[0], model.getParameters()[0]);
}

TEST(model, testScalingLog) {
    model.setParameterScale(ParameterScaling::ln);

    DOUBLES_EQUAL(std::log(p[0]), model.getParameters()[0], 1e-16);
}

TEST(model, testScalingLog10) {
    model.setParameterScale(ParameterScaling::log10);

    DOUBLES_EQUAL(std::log10(p[0]), model.getParameters()[0], 1e-16);
}

TEST(model, testParameterScalingLengthMismatch) {
    // too short
    auto pscale = std::vector<ParameterScaling>(p.size() - 1, ParameterScaling::log10);
    CHECK_THROWS(AmiException, model.setParameterScale(pscale));

    // too long
    pscale = std::vector<ParameterScaling>(p.size() + 1, ParameterScaling::log10);
    CHECK_THROWS(AmiException, model.setParameterScale(pscale));
}

TEST(model, testSetTimepoints) {
    CHECK_THROWS(AmiException,model.setTimepoints(std::vector<realtype>{0.0,1.0,0.5}))
}

TEST(model, testNameIdGetterSetter){
    model.setParameterById("p0",3.0);
    DOUBLES_EQUAL(model.getParameterById("p0"), 3.0, 1e-16);
    CHECK_THROWS(AmiException,model.getParameterById("p1"));
    DOUBLES_EQUAL(model.setParametersByIdRegex("p[\\d]+",5.0), p.size(), 1e-16);
    for (const auto &ip: model.getParameters())
       DOUBLES_EQUAL(ip, 5.0, 1e-16)
    CHECK_THROWS(AmiException,model.setParametersByIdRegex("k[\\d]+", 5.0));
    
    model.setParameterByName("p0",3.0);
    DOUBLES_EQUAL(model.getParameterByName("p0"), 3.0, 1e-16);
    CHECK_THROWS(AmiException,model.getParameterByName("p1"));
    DOUBLES_EQUAL(model.setParametersByNameRegex("p[\\d]+",5.0), p.size(), 1e-16);
    for (const auto &ip: model.getParameters())
        DOUBLES_EQUAL(ip, 5.0, 1e-16)
    CHECK_THROWS(AmiException,model.setParametersByNameRegex("k[\\d]+", 5.0));
    
    model.setFixedParameterById("k0",3.0);
    DOUBLES_EQUAL(model.getFixedParameterById("k0"), 3.0, 1e-16);
    CHECK_THROWS(AmiException,model.getFixedParameterById("k4"));
    DOUBLES_EQUAL(model.setFixedParametersByIdRegex("k[\\d]+",5.0), k.size(), 1e-16);
    for (const auto &ik: model.getFixedParameters())
        DOUBLES_EQUAL(ik, 5.0, 1e-16)
    CHECK_THROWS(AmiException,model.setFixedParametersByIdRegex("p[\\d]+", 5.0));
    
    model.setFixedParameterByName("k0",3.0);
    DOUBLES_EQUAL(model.getFixedParameterByName("k0"),3.0, 1e-16);
    CHECK_THROWS(AmiException,model.getFixedParameterByName("k4"));
    DOUBLES_EQUAL(model.setFixedParametersByNameRegex("k[\\d]+",5.0), k.size(), 1e-16);
    for (const auto &ik: model.getFixedParameters())
        DOUBLES_EQUAL(ik, 5.0, 1e-16)
    CHECK_THROWS(AmiException,model.setFixedParametersByNameRegex("p[\\d]+", 5.0));
}

TEST(model, reinitializeFixedParameterInitialStates){
    CHECK_THROWS(AmiException, model.setReinitializeFixedParameterInitialStates(true));
    model.setReinitializeFixedParameterInitialStates(false);
    CHECK_TRUE(!model.getReinitializeFixedParameterInitialStates());
    AmiVector x(nx);
    AmiVectorArray sx(model.np(),nx);
}

TEST_GROUP(symbolicFunctions)
{
    void setup() {

    }

    void teardown() {

    }
};


TEST(symbolicFunctions, testSign) {
    CHECK_EQUAL(-1, sign(-2));
    CHECK_EQUAL(0, sign(0));
    CHECK_EQUAL(1, sign(2));
}

TEST(symbolicFunctions, testHeaviside) {
    CHECK_EQUAL(0, heaviside(-1));
    CHECK_EQUAL(0, heaviside(0));
    CHECK_EQUAL(1, heaviside(1));
}


TEST(symbolicFunctions, testMin) {
    CHECK_EQUAL(-1, amici::min(-1, 2, 0));
    CHECK_EQUAL(-2, amici::min(1, -2, 0));
    CHECK_TRUE(amici::isNaN(amici::min(amici::getNaN(), amici::getNaN(), 0)));
    CHECK_EQUAL(-1, amici::min(-1, amici::getNaN(), 0));
    CHECK_EQUAL(-1, amici::min(amici::getNaN(), -1, 0));
}

TEST(symbolicFunctions, testMax) {
    CHECK_EQUAL(2, amici::max(-1, 2, 0));
    CHECK_EQUAL(1, amici::max(1, -2, 0));
    CHECK_TRUE(amici::isNaN(amici::max(amici::getNaN(), amici::getNaN(), 0)));
    CHECK_EQUAL(-1, amici::max(-1, amici::getNaN(), 0));
    CHECK_EQUAL(-1, amici::max(amici::getNaN(), -1, 0));
}


TEST(symbolicFunctions, testDMin) {
    CHECK_EQUAL(0, amici::Dmin(1, -1, -2, 0));
    CHECK_EQUAL(1, amici::Dmin(1, -1, 2, 0));
    CHECK_EQUAL(1, amici::Dmin(2, -1, -2, 0));
    CHECK_EQUAL(0, amici::Dmin(2, -1, 2, 0));
}

TEST(symbolicFunctions, testDMax) {
    CHECK_EQUAL(1, amici::Dmax(1, -1, -2, 0));
    CHECK_EQUAL(0, amici::Dmax(1, -1, 2, 0));
    CHECK_EQUAL(0, amici::Dmax(2, -1, -2, 0));
    CHECK_EQUAL(1, amici::Dmax(2, -1, 2, 0));
}

TEST(symbolicFunctions, testpos_pow) {
    CHECK_EQUAL(0, amici::pos_pow(-0.1, 3));
    CHECK_EQUAL(pow(0.1, 3), amici::pos_pow(0.1, 3));
}

TEST_GROUP(amiciSolver)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(amiciSolver, testEquality) {
    IDASolver i1, i2;
    CVodeSolver c1, c2;

    CHECK_TRUE(i1 == i2);
    CHECK_TRUE(c1 == c2);
    CHECK_FALSE(i1 == c1);
}

TEST(amiciSolver, testClone) {
    IDASolver i1;
    auto i2 = std::unique_ptr<Solver>(i1.clone());
    CHECK_TRUE(i1 == *i2);

    CVodeSolver c1;
    auto c2 = std::unique_ptr<Solver>(c1.clone());
    CHECK_TRUE(c1 == *c2);
    CHECK_FALSE(*i2 == *c2);
}


TEST_GROUP(amiciSolverIdas)
{
    void setup() {

    }

    void teardown() {

    }
};

TEST(amiciSolverIdas, testConstructionDestruction) {
    IDASolver solver;
}

TEST_GROUP(edata)
{
    int nx = 1, ny = 2, nz = 3, nmaxevent = 4;
    std::vector<realtype> timepoints = {1, 2, 3, 4};
    
    std::unique_ptr<amici::Model> model = getModel();
    

    Model_Test model_dim = Model_Test(nx, nx, nx, nx, ny, ny, nz, nz,
    nmaxevent, 0, 0, 0, 0, 0, 0, 0, SecondOrderMode::none,
                           std::vector<realtype>(1,0.0),std::vector<realtype>(3,0),std::vector<int>(2,1),
                           std::vector<realtype>(0,0.0),std::vector<int>(0,1));
    void setup() {
        model->setTimepoints(timepoints);
        model->setNMaxEvent(nmaxevent);
        model_dim.setTimepoints(timepoints);
        model_dim.setNMaxEvent(nmaxevent);
    }
    
    void teardown() {
        
    }
};


TEST(edata, testConstructors1) {
    auto edata = ExpData();
    CHECK_TRUE(edata.nytrue() == 0)
    CHECK_TRUE(edata.nztrue() == 0)
    CHECK_TRUE(edata.nmaxevent() == 0)
}
TEST(edata, testConstructors2) {
    auto edata = ExpData(model->nytrue, model->nztrue, model->nMaxEvent());
    CHECK_TRUE(edata.nytrue() == model->nytrue)
    CHECK_TRUE(edata.nztrue() == model->nztrue)
    CHECK_TRUE(edata.nmaxevent() == model->nMaxEvent())
}

TEST(edata, testConstructors3) {
    auto edata = ExpData(model->nytrue, model->nztrue, model->nMaxEvent(), timepoints);
    CHECK_TRUE(edata.nytrue() == model->nytrue)
    CHECK_TRUE(edata.nztrue() == model->nztrue)
    CHECK_TRUE(edata.nmaxevent() == model->nMaxEvent())
    CHECK_TRUE(edata.nt() == model->nt())
    checkEqualArray(timepoints,edata.getTimepoints(), TEST_ATOL, TEST_RTOL, "ts");
}

TEST(edata, testConstructors4) {
    std::vector<realtype> y(ny*timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny*timepoints.size(), 0.1);
    std::vector<realtype> z(nz*nmaxevent, 0.0);
    std::vector<realtype> z_std(nz*nmaxevent, 0.1);
    
    auto edata = ExpData(model_dim.nytrue,
                         model_dim.nztrue,
                         model_dim.nMaxEvent(),
                         timepoints,
                         y,
                         y_std,
                         z,
                         z_std);
    CHECK_TRUE(edata.nytrue() == model_dim.nytrue)
    CHECK_TRUE(edata.nztrue() == model_dim.nztrue)
    CHECK_TRUE(edata.nmaxevent() == model_dim.nMaxEvent())
    CHECK_TRUE(edata.nt() == model_dim.nt())
    checkEqualArray(timepoints,edata.getTimepoints(), TEST_ATOL, TEST_RTOL, "ts");
    checkEqualArray(y,edata.getObservedData(), TEST_ATOL, TEST_RTOL, "observedData");
    checkEqualArray(y_std,edata.getObservedDataStdDev(), TEST_ATOL, TEST_RTOL, "observedDataStdDev");
    checkEqualArray(z,edata.getObservedEvents(), TEST_ATOL, TEST_RTOL, "observedEvents");
    checkEqualArray(z_std,edata.getObservedEventsStdDev(), TEST_ATOL, TEST_RTOL, "observedEventsStdDev");
    
    auto edata_copy = ExpData(edata);
    CHECK_TRUE(edata.nytrue() == edata_copy.nytrue())
    CHECK_TRUE(edata.nztrue() == edata_copy.nztrue())
    CHECK_TRUE(edata.nmaxevent() == edata_copy.nmaxevent())
    CHECK_TRUE(edata.nt() == edata_copy.nt())
    checkEqualArray(edata_copy.getTimepoints(),edata.getTimepoints(), TEST_ATOL, TEST_RTOL, "ts");
    checkEqualArray(edata_copy.getObservedData(),edata.getObservedData(), TEST_ATOL, TEST_RTOL, "observedData");
    checkEqualArray(edata_copy.getObservedDataStdDev(),edata.getObservedDataStdDev(), TEST_ATOL, TEST_RTOL, "observedDataStdDev");
    checkEqualArray(edata_copy.getObservedEvents(),edata.getObservedEvents(), TEST_ATOL, TEST_RTOL, "observedEvents");
    checkEqualArray(edata_copy.getObservedEventsStdDev(),edata.getObservedEventsStdDev(), TEST_ATOL, TEST_RTOL, "observedEventsStdDev");
}

TEST(edata, testConstructors5) {
    model_dim.setTimepoints(timepoints);
    auto edata = ExpData(model_dim);
    CHECK_TRUE(edata.nytrue() == model_dim.nytrue)
    CHECK_TRUE(edata.nztrue() == model_dim.nztrue)
    CHECK_TRUE(edata.nmaxevent() == model_dim.nMaxEvent())
    CHECK_TRUE(edata.nt() == model_dim.nt())
    checkEqualArray(model_dim.getTimepoints(),edata.getTimepoints(), TEST_ATOL, TEST_RTOL, "ts");
}
    
TEST(edata, testDimensionChecks) {

    std::vector<realtype> bad_std(ny, -0.1);
    
    std::vector<realtype> y(ny*timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny*timepoints.size(), 0.1);
    std::vector<realtype> z(nz*nmaxevent, 0.0);
    std::vector<realtype> z_std(nz*nmaxevent, 0.1);
    
    CHECK_THROWS(AmiException,
                 ExpData(model_dim.nytrue,
                         model_dim.nztrue,
                         model_dim.nMaxEvent(),
                         timepoints,
                         z,
                         z_std,
                         z,
                         z_std)
                 )
    
    CHECK_THROWS(AmiException,
                 ExpData(model_dim.nytrue,
                         model_dim.nztrue,
                         model_dim.nMaxEvent(),
                         timepoints,
                         z,
                         bad_std,
                         z,
                         z_std)
                 )
    
    auto edata = ExpData(model_dim);
    
    std::vector<realtype> bad_y(ny*timepoints.size()+1, 0.0);
    std::vector<realtype> bad_y_std(ny*timepoints.size()+1, 0.1);
    std::vector<realtype> bad_z(nz*nmaxevent+1, 0.0);
    std::vector<realtype> bad_z_std(nz*nmaxevent+1, 0.1);
    
    CHECK_THROWS(AmiException,edata.setObservedData(bad_y))
    CHECK_THROWS(AmiException,edata.setObservedDataStdDev(bad_y_std))
    CHECK_THROWS(AmiException,edata.setObservedEvents(bad_z))
    CHECK_THROWS(AmiException,edata.setObservedEventsStdDev(bad_y_std))
    
    std::vector<realtype> bad_single_y(edata.nt()+1, 0.0);
    std::vector<realtype> bad_single_y_std(edata.nt()+1, 0.1);
    std::vector<realtype> bad_single_z(edata.nmaxevent()+1, 0.0);
    std::vector<realtype> bad_single_z_std(edata.nmaxevent()+1, 0.1);
    
    CHECK_THROWS(AmiException,edata.setObservedData(bad_single_y,0))
    CHECK_THROWS(AmiException,edata.setObservedDataStdDev(bad_single_y_std,0))
    CHECK_THROWS(AmiException,edata.setObservedEvents(bad_single_z,0))
    CHECK_THROWS(AmiException,edata.setObservedEventsStdDev(bad_single_y_std,0))
    
    CHECK_THROWS(AmiException,edata.setTimepoints(std::vector<realtype>{0.0,1.0,0.5}))
}

TEST(edata, testSettersGetters) {
    auto edata = ExpData(model_dim);
    
    std::vector<realtype> y(ny*timepoints.size(), 0.0);
    std::vector<realtype> y_std(ny*timepoints.size(), 0.1);
    std::vector<realtype> z(nz*nmaxevent, 0.0);
    std::vector<realtype> z_std(nz*nmaxevent, 0.1);
    
    edata.setObservedData(y);
    checkEqualArray(edata.getObservedData(), y, TEST_ATOL, TEST_RTOL, "ObservedData");
    edata.setObservedDataStdDev(y_std);
    checkEqualArray(edata.getObservedDataStdDev(), y_std, TEST_ATOL, TEST_RTOL, "ObservedDataStdDev");
    edata.setObservedEvents(z);
    checkEqualArray(edata.getObservedEvents(), z, TEST_ATOL, TEST_RTOL, "ObservedEvents");
    edata.setObservedEventsStdDev(z_std);
    checkEqualArray(edata.getObservedEventsStdDev(), z_std, TEST_ATOL, TEST_RTOL, "ObservedEventsStdDev");
    
    std::vector<realtype> single_y(edata.nt(), 0.0);
    std::vector<realtype> single_y_std(edata.nt(), 0.1);
    
    
    for (int iy = 0; iy < ny; ++iy) {
        edata.setObservedData(single_y, iy);
        edata.setObservedDataStdDev(single_y_std, iy);
    }
    CHECK_THROWS(std::exception,edata.setObservedData(single_y, ny))
    CHECK_THROWS(std::exception,edata.setObservedData(single_y, -1))
    CHECK_THROWS(std::exception,edata.setObservedDataStdDev(single_y_std, ny))
    CHECK_THROWS(std::exception,edata.setObservedDataStdDev(single_y_std, -1))
    
    std::vector<realtype> single_z(edata.nmaxevent(), 0.0);
    std::vector<realtype> single_z_std(edata.nmaxevent(), 0.1);
    
    for (int iz = 0; iz < nz; ++iz) {
        edata.setObservedEvents(single_z, iz);
        edata.setObservedEventsStdDev(single_z_std, iz);
    }
    
    CHECK_THROWS(std::exception,edata.setObservedEvents(single_z, nz))
    CHECK_THROWS(std::exception,edata.setObservedEvents(single_z, -1))
    CHECK_THROWS(std::exception,edata.setObservedEventsStdDev(single_z_std, nz))
    CHECK_THROWS(std::exception,edata.setObservedEventsStdDev(single_z_std, -1))
    
    CHECK_TRUE(edata.getObservedDataPtr(0))
    CHECK_TRUE(edata.getObservedDataStdDevPtr(0))
    CHECK_TRUE(edata.getObservedEventsPtr(0))
    CHECK_TRUE(edata.getObservedEventsStdDevPtr(0))
    
    std::vector<realtype> empty(0, 0.0);
    
    edata.setObservedData(empty);
    edata.setObservedDataStdDev(empty);
    edata.setObservedEvents(empty);
    edata.setObservedEventsStdDev(empty);
    
    CHECK_TRUE(!edata.getObservedDataPtr(0))
    CHECK_TRUE(!edata.getObservedDataStdDevPtr(0))
    CHECK_TRUE(!edata.getObservedEventsPtr(0))
    CHECK_TRUE(!edata.getObservedEventsStdDevPtr(0))
    
    checkEqualArray(edata.getObservedData(), empty, TEST_ATOL, TEST_RTOL, "ObservedData");
    checkEqualArray(edata.getObservedDataStdDev(), empty, TEST_ATOL, TEST_RTOL, "ObservedDataStdDev");
    checkEqualArray(edata.getObservedEvents(), empty, TEST_ATOL, TEST_RTOL, "ObservedEvents");
    checkEqualArray(edata.getObservedEventsStdDev(), empty, TEST_ATOL, TEST_RTOL, "ObservedEventsStdDev");
}

TEST_GROUP(solver)
{
    int nx = 1, ny = 2, nz = 3, ne = 0;
    double tol, badtol;
    std::vector<realtype> timepoints = {1, 2, 3, 4};
    
    std::unique_ptr<amici::Model> model = getModel();
    SensitivityMethod sensi_meth;
    SensitivityOrder  sensi;
    int steps, badsteps;
    LinearMultistepMethod lmm;
    NonlinearSolverIteration iter;
    InternalSensitivityMethod ism;
    InterpolationType interp;
    
    
    Model_Test model_dim = Model_Test(nx, nx, nx, nx, ny, ny, nz, nz, ne, 0, 0,
     0, 0, 0, 0, 0, SecondOrderMode::none,
                                      std::vector<realtype>(1,0.0),std::vector<realtype>(3,0),std::vector<int>(2,1),
                                      std::vector<realtype>(0,0.0),std::vector<int>(0,1));
    
    CVodeSolver solver = CVodeSolver();
    
    void setup() {
        tol = 0.01;
        badtol = -0.01;
        sensi_meth = SensitivityMethod::adjoint;
        sensi = SensitivityOrder::first;
        steps = 1000;
        badsteps = -1;
        lmm = LinearMultistepMethod::adams;
        iter = NonlinearSolverIteration::functional;
        ism = InternalSensitivityMethod::staggered1;
        interp = InterpolationType::polynomial;
    }
    
    void teardown() {
        
    }
};

TEST(solver, testSettersGettersNoSetup) {
    testSolverGetterSetters(solver,sensi_meth,sensi,ism,interp,iter,lmm,steps,badsteps,tol,badtol);
}
                
TEST(solver, testSettersGettersWithSetup) {
    
    solver.setSensitivityMethod(sensi_meth);
    CHECK_EQUAL(static_cast<int>(solver.getSensitivityMethod()), static_cast<int>(sensi_meth));
    
    auto rdata = std::unique_ptr<ReturnData>(new ReturnData(solver,&model_dim));
    AmiVector x(nx), dx(nx);
    AmiVectorArray sx(nx,1), sdx(nx,1);
    
    
    model_dim.setInitialStates(std::vector<realtype>{0});
    
    solver.setup(&x,&dx,&sx,&sdx, &model_dim);
    
    testSolverGetterSetters(solver,sensi_meth,sensi,ism,interp,iter,lmm,steps,badsteps,tol,badtol);
}

void testSolverGetterSetters(CVodeSolver solver, SensitivityMethod sensi_meth, SensitivityOrder sensi, InternalSensitivityMethod ism, InterpolationType interp,
                             NonlinearSolverIteration iter, LinearMultistepMethod lmm, int steps, int badsteps, double tol, double badtol) {
    
    solver.setSensitivityMethod(sensi_meth);
    CHECK_EQUAL(static_cast<int>(solver.getSensitivityMethod()), static_cast<int>(sensi_meth));
    
    solver.setSensitivityOrder(sensi);
    CHECK_EQUAL(static_cast<int>(solver.getSensitivityOrder()), static_cast<int>(sensi));
    
    solver.setInternalSensitivityMethod(ism);
    CHECK_EQUAL(static_cast<int>(solver.getInternalSensitivityMethod()), static_cast<int>(ism));
    
    solver.setInterpolationType(interp);
    CHECK_EQUAL(static_cast<int>(solver.getInterpolationType()), static_cast<int>(interp));
    
    solver.setNonlinearSolverIteration(iter);
    CHECK_EQUAL(static_cast<int>(solver.getNonlinearSolverIteration()), static_cast<int>(iter));
    
    solver.setLinearMultistepMethod(lmm);
    CHECK_EQUAL(static_cast<int>(solver.getLinearMultistepMethod ()), static_cast<int>(lmm));
    
    solver.setNewtonPreequilibration(true);
    CHECK_EQUAL(solver.getNewtonPreequilibration(), true);
    
    solver.setStabilityLimitFlag(true);
    CHECK_EQUAL(solver.getStabilityLimitFlag(), true);
    
    CHECK_THROWS(AmiException,solver.setNewtonMaxSteps(badsteps));
    solver.setNewtonMaxSteps(steps);
    CHECK_EQUAL(solver.getNewtonMaxSteps(), steps);
    
    CHECK_THROWS(AmiException,solver.setNewtonMaxLinearSteps(badsteps));
    solver.setNewtonMaxLinearSteps(steps);
    CHECK_EQUAL(solver.getNewtonMaxLinearSteps(), steps);
    
    CHECK_THROWS(AmiException,solver.setMaxSteps(badsteps));
    solver.setMaxSteps(steps);
    CHECK_EQUAL(solver.getMaxSteps(), steps);
    
    CHECK_THROWS(AmiException,solver.setMaxStepsBackwardProblem(badsteps));
    solver.setMaxStepsBackwardProblem(steps);
    CHECK_EQUAL(solver.getMaxStepsBackwardProblem(), steps);
    
    CHECK_THROWS(AmiException,solver.setRelativeTolerance(badtol));
    solver.setRelativeTolerance(tol);
    CHECK_EQUAL(solver.getRelativeTolerance(), tol);
    
    CHECK_THROWS(AmiException,solver.setAbsoluteTolerance(badtol));
    solver.setAbsoluteTolerance(tol);
    CHECK_EQUAL(solver.getAbsoluteTolerance(), tol);
    
    CHECK_THROWS(AmiException,solver.setRelativeToleranceQuadratures(badtol));
    solver.setRelativeToleranceQuadratures(tol);
    CHECK_EQUAL(solver.getRelativeToleranceQuadratures(), tol);
    
    CHECK_THROWS(AmiException,solver.setAbsoluteToleranceQuadratures(badtol));
    solver.setAbsoluteToleranceQuadratures(tol);
    CHECK_EQUAL(solver.getAbsoluteToleranceQuadratures(), tol);
    
    CHECK_THROWS(AmiException,solver.setRelativeToleranceSteadyState(badtol));
    solver.setRelativeToleranceSteadyState(tol);
    CHECK_EQUAL(solver.getRelativeToleranceSteadyState(), tol);
    
    CHECK_THROWS(AmiException,solver.setAbsoluteToleranceSteadyState(badtol));
    solver.setAbsoluteToleranceSteadyState(tol);
    CHECK_EQUAL(solver.getAbsoluteToleranceSteadyState(), tol);
}


