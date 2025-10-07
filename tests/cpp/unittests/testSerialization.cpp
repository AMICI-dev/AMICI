#include "testfunctions.h"
#include <amici/event.h>
#include <amici/model.h>
#include <amici/serialization.h>
#include <amici/solver_cvodes.h>

#include <cmath>
#include <fstream>
#include <gtest/gtest.h>

void checkReturnDataEqual(
    amici::ReturnData const& r, amici::ReturnData const& s
) {
    ASSERT_EQ(r.id, s.id);
    ASSERT_EQ(r.np, s.np);
    ASSERT_EQ(r.nk, s.nk);
    ASSERT_EQ(r.nx_rdata, s.nx_rdata);
    ASSERT_EQ(r.nxtrue_rdata, s.nxtrue_rdata);
    ASSERT_EQ(r.nx_solver, s.nx_solver);
    ASSERT_EQ(r.nx_solver_reinit, s.nx_solver_reinit);
    ASSERT_EQ(r.ny, s.ny);
    ASSERT_EQ(r.nytrue, s.nytrue);
    ASSERT_EQ(r.nz, s.nz);
    ASSERT_EQ(r.nztrue, s.nztrue);
    ASSERT_EQ(r.ne, s.ne);
    ASSERT_EQ(r.nJ, s.nJ);
    ASSERT_EQ(r.nplist, s.nplist);
    ASSERT_EQ(r.nmaxevent, s.nmaxevent);
    ASSERT_EQ(r.nt, s.nt);
    ASSERT_EQ(r.newton_maxsteps, s.newton_maxsteps);
    ASSERT_EQ(r.pscale, s.pscale);
    ASSERT_EQ(static_cast<int>(r.o2mode), static_cast<int>(s.o2mode));
    ASSERT_EQ(static_cast<int>(r.sensi), static_cast<int>(s.sensi));
    ASSERT_EQ(static_cast<int>(r.sensi_meth), static_cast<int>(s.sensi_meth));

    using amici::checkEqualArray;
    checkEqualArray(r.ts, s.ts, 1e-16, 1e-16, "ts");
    checkEqualArray(r.xdot, s.xdot, 1e-16, 1e-16, "xdot");
    checkEqualArray(r.J, s.J, 1e-16, 1e-16, "J");
    checkEqualArray(r.z, s.z, 1e-16, 1e-16, "z");
    checkEqualArray(r.sigmaz, s.sigmaz, 1e-16, 1e-16, "sigmaz");
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

    ASSERT_EQ(r.numsteps, s.numsteps);
    ASSERT_EQ(r.numsteps_b, s.numsteps_b);
    ASSERT_EQ(r.num_rhs_evals, s.num_rhs_evals);
    ASSERT_EQ(r.num_rhs_evals_b, s.num_rhs_evals_b);
    ASSERT_EQ(r.num_err_test_fails, s.num_err_test_fails);
    ASSERT_EQ(r.num_err_test_fails_b, s.num_err_test_fails_b);
    ASSERT_EQ(r.num_non_lin_solv_conv_fails, s.num_non_lin_solv_conv_fails);
    ASSERT_EQ(r.num_non_lin_solv_conv_fails_b, s.num_non_lin_solv_conv_fails_b);
    ASSERT_EQ(r.order, s.order);
    ASSERT_EQ(r.cpu_time, s.cpu_time);
    ASSERT_EQ(r.cpu_time_b, s.cpu_time_b);

    ASSERT_EQ(r.preeq_status, s.preeq_status);
    ASSERT_TRUE(
        r.preeq_t == s.preeq_t
        || (std::isnan(r.preeq_t) && std::isnan(s.preeq_t))
    );
    ASSERT_TRUE(
        r.preeq_wrms == s.preeq_wrms
        || (std::isnan(r.preeq_wrms) && std::isnan(s.preeq_wrms))
    );
    ASSERT_EQ(r.preeq_numsteps, s.preeq_numsteps);
    EXPECT_NEAR(r.preeq_cpu_time, s.preeq_cpu_time, 1e-16);

    ASSERT_EQ(r.posteq_status, s.posteq_status);
    ASSERT_TRUE(
        r.posteq_t == s.posteq_t
        || (std::isnan(r.posteq_t) && std::isnan(s.posteq_t))
    );
    ASSERT_TRUE(
        r.posteq_wrms == s.posteq_wrms
        || (std::isnan(r.posteq_wrms) && std::isnan(s.posteq_wrms))
    );
    ASSERT_EQ(r.posteq_numsteps, s.posteq_numsteps);
    EXPECT_NEAR(r.posteq_cpu_time, s.posteq_cpu_time, 1e-16);

    checkEqualArray(r.x0, s.x0, 1e-16, 1e-16, "x0");
    checkEqualArray(r.sx0, s.sx0, 1e-16, 1e-16, "sx0");

    ASSERT_TRUE(r.llh == s.llh || (std::isnan(r.llh) && std::isnan(s.llh)));
    ASSERT_TRUE(r.chi2 == s.chi2 || (std::isnan(r.llh) && std::isnan(s.llh)));
    ASSERT_EQ(r.status, s.status);

    checkEqualArray(r.sllh, s.sllh, 1e-5, 1e-5, "sllh");
    checkEqualArray(r.s2llh, s.s2llh, 1e-5, 1e-5, "s2llh");
}

class SolverSerializationTest : public ::testing::Test {
  protected:
    void SetUp() override {
        // set non-default values for all members
        solver.set_absolute_tolerance(1e-4);
        solver.set_relative_tolerance(1e-5);
        solver.set_absolute_tolerance_quadratures(1e-6);
        solver.set_relative_tolerance_quadratures(1e-7);
        solver.set_absolute_tolerance_steady_state(1e-8);
        solver.set_relative_tolerance_steady_state(1e-9);
        solver.set_sensitivity_method(amici::SensitivityMethod::adjoint);
        solver.set_sensitivity_order(amici::SensitivityOrder::second);
        solver.set_max_steps(1e1);
        solver.set_max_steps_backward_problem(1e2);
        solver.set_newton_max_steps(1e3);
        solver.set_state_ordering(
            static_cast<int>(amici::SUNLinSolKLU::StateOrdering::COLAMD)
        );
        solver.set_interpolation_type(amici::InterpolationType::polynomial);
        solver.set_stability_limit_flag(false);
        solver.set_linear_solver(amici::LinearSolver::dense);
        solver.set_linear_multistep_method(amici::LinearMultistepMethod::adams);
        solver.set_non_linear_solver_iteration(
            amici::NonlinearSolverIteration::newton
        );
        solver.set_internal_sensitivity_method(
            amici::InternalSensitivityMethod::staggered
        );
        solver.set_return_data_reporting_mode(
            amici::RDataReporting::likelihood
        );
    }

    amici::CVodeSolver solver;
};

TEST(ModelSerializationTest, ToFile) {
    using amici::realtype;
    int np = 1;
    int nk = 2;
    int nx = 3;
    int ny = 4;
    int nz = 5;
    int ne = 6;
    amici::CVodeSolver solver;
    amici::Model_Test m = amici::Model_Test(
        amici::ModelDimensions{
            .nx_rdata = nx,
            .nxtrue_rdata = nx,
            .nx_solver = nx,
            .nxtrue_solver = nx,
            .nx_solver_reinit = 0,
            .np = np,
            .nk = nk,
            .ny = ny,
            .nytrue = ny,
            .nz = nz,
            .nztrue = nz,
            .ne = ne,
            .ne_solver = 0,
            .nspl = 0,
            .nw = 9,
            .ndwdx = 2,
            .ndwdp = 2,
            .ndwdw = 2,
            .ndxdotdw = 13,
            .ndJydy = {0, 0, 0, 0},
            .ndxrdatadxsolver = 9,
            .ndxrdatadtcl = 0,
            .ndtotal_cldx_rdata = 0,
            .nnz = 17,
            .nJ = 0,
            .ubw = 18,
            .lbw = 19
        },

        amici::SimulationParameters(
            std::vector<realtype>(nk, 0.0), std::vector<realtype>(np, 0.0),
            std::vector<int>(np, 0)
        ),
        amici::SecondOrderMode::none, std::vector<realtype>(nx, 0.0),
        std::vector<int>(nz, 0),
        {
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
        }
    );

    {
        std::ofstream ofs("sstore.dat");
        boost::archive::text_oarchive oar(ofs);
        // oar & static_cast<amici::Solver&>(solver);
        oar& static_cast<amici::Model&>(m);
    }
    {
        std::ifstream ifs("sstore.dat");
        boost::archive::text_iarchive iar(ifs);
        amici::CVodeSolver v;
        amici::Model_Test n;
        // iar &static_cast<amici::Solver&>(v);
        iar& static_cast<amici::Model&>(n);
        // CHECK_TRUE(solver == v);
        ASSERT_EQ(m, n);
    }
}

TEST(ReturnDataSerializationTest, ToString) {
    using amici::realtype;
    int np = 1;
    int nk = 2;
    int nx = 3;
    int ny = 4;
    int nz = 5;
    int ne = 6;
    amici::CVodeSolver solver;
    amici::Model_Test m = amici::Model_Test(
        amici::ModelDimensions{
            .nx_rdata = nx,
            .nxtrue_rdata = nx,
            .nx_solver = nx,
            .nxtrue_solver = nx,
            .nx_solver_reinit = 0,
            .np = np,
            .nk = nk,
            .ny = ny,
            .nytrue = ny,
            .nz = nz,
            .nztrue = nz,
            .ne = ne,
            .ne_solver = 0,
            .nspl = 0,
            .nw = 9,
            .ndwdx = 10,
            .ndwdp = 2,
            .ndwdw = 12,
            .ndxdotdw = 13,
            .ndJydy = {0, 0, 0, 0},
            .ndxrdatadxsolver = 9,
            .ndxrdatadtcl = 0,
            .ndtotal_cldx_rdata = 0,
            .nnz = 17,
            .nJ = 0,
            .ubw = 18,
            .lbw = 19
        },

        amici::SimulationParameters(
            std::vector<realtype>(nk, 0.0), std::vector<realtype>(np, 0.0),
            std::vector<int>(np, 0)
        ),
        amici::SecondOrderMode::none, std::vector<realtype>(nx, 0.0),
        std::vector<int>(nz, 0),
        {
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
            amici::Event("e1", true, true, 1),
        }
    );

    amici::ReturnData r(solver, m);
    r.id = "some_id";
    std::string serialized = amici::serialize_to_string(r);

    checkReturnDataEqual(
        r, amici::deserialize_from_string<amici::ReturnData>(serialized)
    );
}

TEST_F(SolverSerializationTest, ToChar) {
    int length;
    char* buf = amici::serialize_to_char(solver, &length);

    amici::CVodeSolver v
        = amici::deserialize_from_char<amici::CVodeSolver>(buf, length);

    delete[] buf;
    ASSERT_EQ(solver, v);
}

TEST_F(SolverSerializationTest, ToStdVec) {

    auto buf = amici::serialize_to_std_vec(solver);
    amici::CVodeSolver v = amici::deserialize_from_char<amici::CVodeSolver>(
        buf.data(), buf.size()
    );

    ASSERT_EQ(solver, v);
}
