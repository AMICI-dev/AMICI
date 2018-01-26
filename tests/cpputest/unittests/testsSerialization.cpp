
#include <include/amici_serialization.h>
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <include/amici_model.h>

void checkUserDataEqual(amici::UserData const& u, amici::UserData const& v) {
    CHECK_EQUAL(u.np(), v.np());
    CHECK_EQUAL(u.nx(), v.nx());
    CHECK_EQUAL(u.nk(), v.nk());
    //CHECK_EQUAL(u.pscale, v.pscale);

    CHECK_EQUAL(u.nme(), v.nme());
    CHECK_EQUAL(u.nplist(), v.nplist());
    CHECK_EQUAL(u.nt(), v.nt());
    CHECK_EQUAL(u.t0(), v.t0());
    CHECK_EQUAL(u.sensi, v.sensi);
    CHECK_EQUAL(u.atol, v.atol);
    CHECK_EQUAL(u.rtol, v.rtol);
    CHECK_EQUAL(u.maxsteps, v.maxsteps);
    CHECK_EQUAL(u.ism, v.ism);
    CHECK_EQUAL(u.sensmeth(), v.sensmeth());
    CHECK_EQUAL(u.linsol, v.linsol);
    //CHECK_EQUAL(u.interpType, v.interpType);
    //CHECK_EQUAL(u.lmm, v.lmm);
    //CHECK_EQUAL(u.iter, v.iter);
    //CHECK_EQUAL(u.stldet, v.stldet);
    //CHECK_EQUAL(u.ordering, v.ordering);

    amici::checkEqualArray(u.p(), v.p(), u.np(), 1e-16, 1e-16, "p");
    amici::checkEqualArray(u.k(), v.k(), u.nk(), 1e-16, 1e-16, "k");
}


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

    checkEqualArray(r.sllh, s.sllh, r.nplist, 1e-8, 1e-8, "sllh");
    checkEqualArray(r.s2llh, s.s2llh, r.nplist * r.nplist, 1e-8, 1e-8, "s2llh");
}


TEST_GROUP(userDataSerialization){void setup(){}

                          void teardown(){}};


TEST(userDataSerialization, testFile) {

    amici::UserData u(1, 2, 3);
    {
        std::ofstream ofs("sstore.dat");
        boost::archive::text_oarchive oar(ofs);
        oar &u;
    }
    {
        std::ifstream ifs("sstore.dat");
        boost::archive::text_iarchive iar(ifs);
        amici::UserData v;
        iar &v;
        checkUserDataEqual(u, v);
    }
}

TEST(userDataSerialization, testString) {
    amici::UserData u(1, 2, 3);

    std::string serialized = serializeToString(u);

    checkUserDataEqual(u, amici::deserializeFromString<amici::UserData>(serialized));
}

TEST(userDataSerialization, testChar) {
    amici::UserData u(2, 1, 3);
    double p[2] = {1,2};
    u.setParameters(p);

    int length;
    char *buf = serializeToChar(&u, &length);

    amici::UserData v = amici::deserializeFromChar<amici::UserData>(buf, length);

    delete[] buf;
    checkUserDataEqual(u, v);
}

TEST(userDataSerialization, testStdVec) {

    amici::UserData u(2, 1, 3);
    double p[2] = {1,2};
    u.setParameters(p);

    auto buf = amici::serializeToStdVec(&u);

    amici::UserData v = amici::deserializeFromChar<amici::UserData>(buf.data(), buf.size());

    checkUserDataEqual(u, v);
}


TEST_GROUP(returnDataSerialization){void setup(){}

                          void teardown(){}};

TEST(returnDataSerialization, testString) {
    int np = 1;
    int nk = 2;
    int nx = 3;
    int nz = 4;
    amici::UserData u(np, nx, nk);
    amici::Model_Test m( nx, nx, 4, 4, nz, nz, 8, 9, 10, 11, 12, 13, 14, 15, amici::AMICI_O2MODE_NONE,
                   std::vector<realtype>(np,0.0),std::vector<realtype>(nk,0.0),std::vector<int>(np,0),
                   std::vector<realtype>(nx,0.0),std::vector<int>(nz,0));

    amici::ReturnData r(&u, &m);

    std::string serialized = amici::serializeToString(r);

    checkReturnDataEqual(r, amici::deserializeFromString<amici::ReturnData>(serialized));
}

