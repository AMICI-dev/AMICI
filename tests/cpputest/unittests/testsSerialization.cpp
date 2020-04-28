#include <amici/serialization.h> // needs to be included before cpputest
#include <amici/model.h>
#include <amici/solver_cvodes.h>

#include "testfunctions.hpp"

#include <cmath>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

void checkReturnDataEqual(amici::ReturnData const& r, amici::ReturnData const& s) {
    CHECK_EQUAL(r.np, s.np);
    CHECK_EQUAL(r.nk, s.nk);
    CHECK_EQUAL(r.nx, s.nx);
    CHECK_EQUAL(r.nxtrue, s.nxtrue);
    CHECK_EQUAL(r.nx_solver, s.nx_solver);
    CHECK_EQUAL(r.ny, s.ny);
    CHECK_EQUAL(r.nytrue, s.nytrue);
    CHECK_EQUAL(r.nz, s.nz);
    CHECK_EQUAL(r.nztrue, s.nztrue);
    CHECK_EQUAL(r.ne, s.ne);
    CHECK_EQUAL(r.nJ, s.nJ);
    CHECK_EQUAL(r.nplist, s.nplist);
    CHECK_EQUAL(r.nmaxevent, s.nmaxevent);
    CHECK_EQUAL(r.nt, s.nt);
    CHECK_EQUAL(r.newton_maxsteps, s.newton_maxsteps);
    CHECK_TRUE(r.pscale == s.pscale);
    CHECK_EQUAL(static_cast<int>(r.o2mode), static_cast<int>(s.o2mode));
    CHECK_EQUAL(static_cast<int>(r.sensi), static_cast<int>(s.sensi));
    CHECK_EQUAL(static_cast<int>(r.sensi_meth), static_cast<int>(s.sensi_meth));

    using amici::checkEqualArray;
    checkEqualArray(r.ts, s.ts, 1e-16, 1e-16, "ts");
    checkEqualArray(r.xdot, s.xdot, 1e-16, 1e-16, "xdot");
    checkEqualArray(r.J, s.J, 1e-16, 1e-16, "J");
    checkEqualArray(r.z, s.z, 1e-16, 1e-16, "z");
    checkEqualArray(r.sigmaz, s.sigmaz,1e-16, 1e-16, "sigmaz");
    checkEqualArray(r.sz, s.sz, 1e-16, 1e-16, "sz");
    checkEqualArray(r.ssigmaz, s.ssigmaz, 1e-16, 1e-16, "ssigmaz");
    checkEqualArray(r.rz, s.rz, 1e-16, 1e-16, "rz");
    checkEqualArray(r.srz, s.srz, 1e-16, 1e-16, "srz");
    checkEqualArray(r.s2rz, s.s2rz, 1e-16, 1e-16, "s2rz");
    checkEqualArray(r.x, s.x, 1e-16, 1e-16, "x");
    checkEqualArray(r.sx, s.sx, 1e-16, 1e-16, "sx");

    checkEqualArray(r.y, s.y, 1e-16, 1e-16, "y");
    checkEqualArray(r.sigmay, s.sigmay, 1e-16, 1e-16, "sigmay");
    checkEqualArray(r.sy, s.sy, 1e-16, 1e-16, "sy");
    checkEqualArray(r.ssigmay, s.ssigmay, 1e-16, 1e-16, "ssigmay");

    CHECK_TRUE(r.numsteps == s.numsteps);
    CHECK_TRUE(r.numstepsB == s.numstepsB);
    CHECK_TRUE(r.numrhsevals == s.numrhsevals);
    CHECK_TRUE(r.numrhsevalsB == s.numrhsevalsB);
    CHECK_TRUE(r.numerrtestfails == s.numerrtestfails);
    CHECK_TRUE(r.numerrtestfailsB == s.numerrtestfailsB);
    CHECK_TRUE(r.numnonlinsolvconvfails == s.numnonlinsolvconvfails);
    CHECK_TRUE(r.numnonlinsolvconvfailsB == s.numnonlinsolvconvfailsB);
    CHECK_TRUE(r.order == s.order);
    CHECK_TRUE(r.cpu_time == s.cpu_time);
    CHECK_TRUE(r.cpu_timeB == s.cpu_timeB);

    CHECK_TRUE(r.preeq_status == s.preeq_status);
    CHECK_TRUE(r.preeq_t == s.preeq_t ||
               (std::isnan(r.preeq_t) && std::isnan(s.preeq_t)));
    CHECK_TRUE(r.preeq_wrms == s.preeq_wrms ||
               (std::isnan(r.preeq_wrms) && std::isnan(s.preeq_wrms)));
    CHECK_TRUE(r.preeq_numsteps == s.preeq_numsteps);
    CHECK_TRUE(r.preeq_numlinsteps == s.preeq_numlinsteps);
    DOUBLES_EQUAL(r.preeq_cpu_time, s.preeq_cpu_time, 1e-16);

    CHECK_TRUE(r.posteq_status == s.posteq_status);
    CHECK_TRUE(r.posteq_t == s.posteq_t ||
               (std::isnan(r.posteq_t) && std::isnan(s.posteq_t)));
    CHECK_TRUE(r.posteq_wrms == s.posteq_wrms ||
               (std::isnan(r.posteq_wrms) && std::isnan(s.posteq_wrms)));
    CHECK_TRUE(r.posteq_numsteps == s.posteq_numsteps);
    CHECK_TRUE(r.posteq_numlinsteps == s.posteq_numlinsteps);
    DOUBLES_EQUAL(r.posteq_cpu_time, s.posteq_cpu_time, 1e-16);

    checkEqualArray(r.x0, s.x0, 1e-16, 1e-16, "x0");
    checkEqualArray(r.sx0, s.sx0, 1e-16, 1e-16, "sx0");

    CHECK_TRUE(r.llh == s.llh || (std::isnan(r.llh) && std::isnan(s.llh)));
    CHECK_TRUE(r.chi2 == s.chi2 || (std::isnan(r.llh) && std::isnan(s.llh)));
    CHECK_EQUAL(r.status, s.status);

    checkEqualArray(r.sllh, s.sllh, 1e-5, 1e-5, "sllh");
    checkEqualArray(r.s2llh, s.s2llh, 1e-5, 1e-5, "s2llh");
}


// clang-format off
TEST_GROUP(dataSerialization){
    amici::CVodeSolver solver;
    void setup() {
        // set non-default values for all members
        solver.setAbsoluteTolerance(1e-4);
        solver.setRelativeTolerance(1e-5);
        solver.setAbsoluteToleranceQuadratures(1e-6);
        solver.setRelativeToleranceQuadratures(1e-7);
        solver.setAbsoluteToleranceSteadyState(1e-8);
        solver.setRelativeToleranceSteadyState(1e-9);
        solver.setSensitivityMethod(amici::SensitivityMethod::adjoint);
        solver.setSensitivityOrder(amici::SensitivityOrder::second);
        solver.setMaxSteps(1e1);
        solver.setMaxStepsBackwardProblem(1e2);
        solver.setNewtonMaxSteps(1e3);
        solver.setNewtonMaxLinearSteps(1e4);
        solver.setPreequilibration(true);
        solver.setStateOrdering(static_cast<int>(amici::SUNLinSolKLU::StateOrdering::COLAMD));
        solver.setInterpolationType(amici::InterpolationType::polynomial);
        solver.setStabilityLimitFlag(0);
        solver.setLinearSolver(amici::LinearSolver::dense);
        solver.setLinearMultistepMethod(amici::LinearMultistepMethod::adams);
        solver.setNonlinearSolverIteration(amici::NonlinearSolverIteration::newton);
        solver.setInternalSensitivityMethod(amici::InternalSensitivityMethod::staggered);
        solver.setReturnDataReportingMode(amici::RDataReporting::likelihood);
    }

    void teardown() {
    }
};
// clang-format on



TEST(dataSerialization, testFile) {
    int np = 1;
    int nk = 2;
    int nx = 3;
    int nz = 4;
    amici::CVodeSolver solver;
    amici::Model_Test m = amici::Model_Test(nx, nx, nx, nx, 4, 4, nz, nz, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                                            amici::SecondOrderMode::none,
                                            std::vector<realtype>(np,0.0),
                                            std::vector<realtype>(nk,0.0),
                                            std::vector<int>(np,0),
                                            std::vector<realtype>(nx,0.0),
                                            std::vector<int>(nz,0));

    {
        std::ofstream ofs("sstore.dat");
        boost::archive::text_oarchive oar(ofs);
        //oar & static_cast<amici::Solver&>(solver);
        oar & static_cast<amici::Model&>(m);
    }
    {
        std::ifstream ifs("sstore.dat");
        boost::archive::text_iarchive iar(ifs);
        amici::CVodeSolver v;
        amici::Model_Test n;
        //iar &static_cast<amici::Solver&>(v);
        iar &static_cast<amici::Model&>(n);
        //CHECK_TRUE(solver == v);
        CHECK_TRUE(m == n);

    }
}

TEST(dataSerialization, testString) {
    int np = 1;
    int nk = 2;
    int nx = 3;
    int nz = 4;
    amici::CVodeSolver solver;
    amici::Model_Test m = amici::Model_Test(nx, nx, nx, nx, 4, 4, nz, nz, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                                            amici::SecondOrderMode::none,
                                            std::vector<realtype>(np,0.0),
                                            std::vector<realtype>(nk,0.0),
                                            std::vector<int>(np,0),
                                            std::vector<realtype>(nx,0.0),
                                            std::vector<int>(nz,0));

    amici::ReturnData r(solver, m);

    std::string serialized = amici::serializeToString(r);

    checkReturnDataEqual(r, amici::deserializeFromString<amici::ReturnData>(serialized));
}

TEST(dataSerialization, testChar) {
    int length;
    char *buf = amici::serializeToChar(solver, &length);

    amici::CVodeSolver v = amici::deserializeFromChar<amici::CVodeSolver>(buf, length);

    delete[] buf;
    CHECK_TRUE(solver == v);
}

TEST(dataSerialization, testStdVec) {

    auto buf = amici::serializeToStdVec(solver);
    amici::CVodeSolver v = amici::deserializeFromChar<amici::CVodeSolver>(buf.data(), buf.size());

    CHECK_TRUE(solver == v);
}
