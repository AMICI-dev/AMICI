#include <include/amici_serialization.h>
#include "testfunctions.h"
#include <include/amici_model.h>
#include <include/amici_solver_cvodes.h>

#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"


void checkReturnDataEqual(amici::ReturnData const& r, amici::ReturnData const& s) {
    CHECK_EQUAL(r.np, s.np);
    CHECK_EQUAL(r.nk, s.nk);
    CHECK_EQUAL(r.nx, s.nx);
    CHECK_EQUAL(r.nxtrue, s.nxtrue);
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
    CHECK_EQUAL(r.pscale, s.pscale);
    CHECK_EQUAL(r.o2mode, s.o2mode);
    CHECK_EQUAL(r.sensi, s.sensi);
    CHECK_EQUAL(r.sensi_meth, s.sensi_meth);

    using amici::checkEqualArray;
    checkEqualArray(r.ts, s.ts, r.nt, 1e-16, 1e-16, "ts");
    checkEqualArray(r.xdot, s.xdot, r.nx, 1e-16, 1e-16, "xdot");
    checkEqualArray(r.J, s.J, r.nx * r.nx, 1e-16, 1e-16, "J");
    checkEqualArray(r.z, s.z, r.nmaxevent * r.nz, 1e-16, 1e-16, "z");
    checkEqualArray(r.sigmaz, s.sigmaz, r.nmaxevent * r.nz, 1e-16, 1e-16, "sigmaz");
    checkEqualArray(r.sz, s.sz, r.nmaxevent * r.nz * r.nplist, 1e-16, 1e-16, "sz");
    checkEqualArray(r.ssigmaz, s.ssigmaz, r.nmaxevent * r.nz * r.nplist, 1e-16, 1e-16, "ssigmaz");
    checkEqualArray(r.rz, s.rz, r.nmaxevent * r.nz, 1e-16, 1e-16, "rz");
    checkEqualArray(r.srz, s.srz, r.nmaxevent * r.nz * r.nplist, 1e-16, 1e-16, "srz");
    checkEqualArray(r.s2rz, s.s2rz, r.nmaxevent * r.nz * r.nplist * r.nplist, 1e-16, 1e-16, "s2rz");
    checkEqualArray(r.x, s.x, r.nt * r.nx, 1e-16, 1e-16, "x");
    checkEqualArray(r.sx, s.sx, r.nt * r.nx * r.nplist, 1e-16, 1e-16, "sx");

    checkEqualArray(r.y, s.y, r.nt * r.ny, 1e-16, 1e-16, "y");
    checkEqualArray(r.sigmay, s.sigmay, r.nt * r.ny * r.nplist, 1e-16, 1e-16, "sigmay");
    checkEqualArray(r.sy, s.sy, r.nt * r.ny * r.nplist, 1e-16, 1e-16, "sy");
    checkEqualArray(r.ssigmay, s.ssigmay, r.nt * r.ny * r.nplist, 1e-16, 1e-16, "ssigmay");

    checkEqualArray(r.numsteps, s.numsteps, r.nt, 1e-16, 1e-16, "numsteps");
    checkEqualArray(r.numstepsB, s.numstepsB, r.nt, 1e-16, 1e-16, "numstepsB");
    checkEqualArray(r.numrhsevals, s.numrhsevals, r.nt, 1e-16, 1e-16, "numrhsevals");
    checkEqualArray(r.numrhsevalsB, s.numrhsevalsB, r.nt, 1e-16, 1e-16, "numrhsevalsB");
    checkEqualArray(r.numerrtestfails, s.numerrtestfails, r.nt, 1e-16, 1e-16, "numerrtestfails");
    checkEqualArray(r.numerrtestfailsB, s.numerrtestfailsB, r.nt, 1e-16, 1e-16, "numerrtestfailsB");
    checkEqualArray(r.numnonlinsolvconvfails, s.numnonlinsolvconvfails, r.nt, 1e-16, 1e-16, "numnonlinsolvconvfails");
    checkEqualArray(r.numnonlinsolvconvfailsB, s.numnonlinsolvconvfailsB, r.nt, 1e-16, 1e-16, "numnonlinsolvconvfailsB");
    checkEqualArray(r.order, s.order, r.nt, 1e-16, 1e-16, "order");

    checkEqualArray(r.newton_status, s.newton_status, r.nt, 1e-16, 1e-16, "newton_status");
    checkEqualArray(r.newton_time, s.newton_time, r.nt, 1e-16, 1e-16, "newton_time");
    checkEqualArray(r.newton_numsteps, s.newton_numsteps, r.nt, 1e-16, 1e-16, "newton_numsteps");
    checkEqualArray(r.newton_numlinsteps, s.newton_numlinsteps, r.nt, 1e-16, 1e-16, "newton_numlinsteps");
    checkEqualArray(r.x0, s.x0, r.nx, 1e-16, 1e-16, "x0");
    checkEqualArray(r.sx0, s.sx0, r.nx * r.nplist, 1e-16, 1e-16, "sx0");

    CHECK_EQUAL(*r.llh, *s.llh);
    CHECK_EQUAL(*r.chi2, *s.chi2);
    CHECK_EQUAL(*r.status, *s.status);

    checkEqualArray(r.sllh, s.sllh, r.nplist, 1e-5, 1e-5, "sllh");
    checkEqualArray(r.s2llh, s.s2llh, r.nplist * r.nplist, 1e-5, 1e-5, "s2llh");
}


// clang-format off
TEST_GROUP(dataSerialization){

    void setup() {

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
    amici::Model_Test m = amici::Model_Test(nx, nx, 4, 4, nz, nz, 8, 9, 10, 11, 12, 13, 14, 15, amici::AMICI_O2MODE_NONE,
                                             std::vector<realtype>(np,0.0),std::vector<realtype>(nk,0.0),std::vector<int>(np,0),
                                             std::vector<realtype>(nx,0.0),std::vector<int>(nz,0));

    {
        std::ofstream ofs("sstore.dat");
        boost::archive::text_oarchive oar(ofs);
        oar & static_cast<amici::Solver&>(solver);
        oar & static_cast<amici::Model&>(m);
    }
    {
        std::ifstream ifs("sstore.dat");
        boost::archive::text_iarchive iar(ifs);
        amici::CVodeSolver v;
        amici::Model_Test n;
        iar &static_cast<amici::Solver&>(v);
        iar &static_cast<amici::Model&>(n);
        CHECK_TRUE(solver == v);
        CHECK_TRUE(m == n);

    }
}

TEST(dataSerialization, testChar) {
    amici::CVodeSolver solver;
    solver.setAbsoluteTolerance(4);

    int length;
    char *buf = amici::serializeToChar<amici::CVodeSolver>(&solver, &length);

    amici::CVodeSolver v = amici::deserializeFromChar<amici::CVodeSolver>(buf, length);

    delete[] buf;
    CHECK_TRUE(solver == v);
}

TEST(dataSerialization, testStdVec) {
    amici::CVodeSolver solver;
    solver.setAbsoluteTolerance(4);

    auto buf = amici::serializeToStdVec<amici::CVodeSolver>(&solver);
    amici::CVodeSolver v = amici::deserializeFromChar<amici::CVodeSolver>(buf.data(), buf.size());

    CHECK_TRUE(solver == v);
}


TEST(dataSerialization, testString) {
    int np = 1;
    int nk = 2;
    int nx = 3;
    int nz = 4;
    amici::CVodeSolver solver;
    amici::Model_Test m = amici::Model_Test(nx, nx, 4, 4, nz, nz, 8, 9, 10, 11, 12, 13, 14, 15, amici::AMICI_O2MODE_NONE,
                                             std::vector<realtype>(np,0.0),std::vector<realtype>(nk,0.0),std::vector<int>(np,0),
                                             std::vector<realtype>(nx,0.0),std::vector<int>(nz,0));

    amici::ReturnData r(solver, &m);

    std::string serialized = amici::serializeToString(r);

    checkReturnDataEqual(r, amici::deserializeFromString<amici::ReturnData>(serialized));
}

