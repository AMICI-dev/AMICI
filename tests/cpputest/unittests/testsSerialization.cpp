
#include <include/amici_serialization.h>
#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include "testfunctions.h"

#include <include/amici_model.h>

void checkUserDataEqual(UserData const& u, UserData const& v) {
    CHECK_EQUAL(u.np, v.np);
    CHECK_EQUAL(u.nx, v.nx);
    CHECK_EQUAL(u.nk, v.nk);
    CHECK_EQUAL(u.nx, v.nx);
    CHECK_EQUAL(u.pscale, v.pscale);

    CHECK_EQUAL(u.nmaxevent, v.nmaxevent);
    CHECK_EQUAL(u.nplist, v.nplist);
    CHECK_EQUAL(u.nt, v.nt);
    CHECK_EQUAL(u.tstart, v.tstart);
    CHECK_EQUAL(u.sensi, v.sensi);
    CHECK_EQUAL(u.atol, v.atol);
    CHECK_EQUAL(u.rtol, v.rtol);
    CHECK_EQUAL(u.maxsteps, v.maxsteps);
    CHECK_EQUAL(u.ism, v.ism);
    CHECK_EQUAL(u.sensi_meth, v.sensi_meth);
    CHECK_EQUAL(u.linsol, v.linsol);
    CHECK_EQUAL(u.interpType, v.interpType);
    CHECK_EQUAL(u.lmm, v.lmm);
    CHECK_EQUAL(u.iter, v.iter);
    CHECK_EQUAL(u.stldet, v.stldet);
    CHECK_EQUAL(u.ordering, v.ordering);

    checkEqualArray(u.p, v.p, u.np, 1e-16, 1e-16);
    checkEqualArray(u.k, v.k, u.nk, 1e-16, 1e-16);
}


void checkReturnDataEqual(ReturnData const& r, ReturnData const& s) {
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

    checkEqualArray(r.ts, s.ts, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.xdot, s.xdot, r.nx, 1e-16, 1e-16);
    checkEqualArray(r.J, s.J, r.nx * r.nx, 1e-16, 1e-16);
    checkEqualArray(r.z, s.z, r.nmaxevent * r.nz, 1e-16, 1e-16);
    checkEqualArray(r.sigmaz, s.sigmaz, r.nmaxevent * r.nz, 1e-16, 1e-16);
    checkEqualArray(r.sz, s.sz, r.nmaxevent * r.nz * r.nplist, 1e-16, 1e-16);
    checkEqualArray(r.ssigmaz, s.sigmaz, r.nmaxevent * r.nz * r.nplist, 1e-16, 1e-16);
    checkEqualArray(r.rz, s.rz, r.nmaxevent * r.nz, 1e-16, 1e-16);
    checkEqualArray(r.srz, s.srz, r.nmaxevent * r.nz * r.nplist, 1e-16, 1e-16);
    checkEqualArray(r.s2rz, s.s2rz, r.nmaxevent * r.nz * r.nplist * r.nplist, 1e-16, 1e-16);
    checkEqualArray(r.ssigmaz, s.sigmaz, r.nmaxevent * r.nz * r.nplist, 1e-16, 1e-16);
    checkEqualArray(r.x, s.x, r.nt * r.nx, 1e-16, 1e-16);
    checkEqualArray(r.sx, s.sx, r.nt * r.nx * r.nplist, 1e-16, 1e-16);

    checkEqualArray(r.y, s.y, r.nt * r.ny, 1e-16, 1e-16);
    checkEqualArray(r.sigmay, s.sigmay, r.nt * r.ny * r.nplist, 1e-16, 1e-16);
    checkEqualArray(r.sy, s.sy, r.nt * r.ny * r.nplist, 1e-16, 1e-16);
    checkEqualArray(r.ssigmay, s.ssigmay, r.nt * r.ny * r.nplist, 1e-16, 1e-16);

    checkEqualArray(r.numsteps, s.numsteps, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.numstepsB, s.numstepsB, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.numrhsevals, s.numrhsevals, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.numrhsevalsB, s.numrhsevalsB, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.numerrtestfails, s.numerrtestfails, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.numerrtestfailsB, s.numerrtestfailsB, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.numnonlinsolvconvfails, s.numnonlinsolvconvfails, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.numnonlinsolvconvfailsB, s.numnonlinsolvconvfailsB, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.order, s.order, r.nt, 1e-16, 1e-16);

    checkEqualArray(r.newton_status, s.newton_status, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.newton_time, s.newton_time, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.newton_numsteps, s.newton_numsteps, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.newton_numlinsteps, s.newton_numlinsteps, r.nt, 1e-16, 1e-16);
    checkEqualArray(r.x0, s.x0, r.nx, 1e-16, 1e-16);
    checkEqualArray(r.sx0, s.sx0, r.nx * r.nplist, 1e-16, 1e-16);

    CHECK_EQUAL(*r.llh, *s.llh);
    CHECK_EQUAL(*r.chi2, *s.chi2);
    CHECK_EQUAL(*r.status, *s.status);

    checkEqualArray(r.sllh, s.sllh, r.nplist, 1e-16, 1e-16);
    checkEqualArray(r.s2llh, s.s2llh, r.nplist * r.nplist, 1e-16, 1e-16);
}


TEST_GROUP(userDataSerialization){void setup(){}

                          void teardown(){}};


TEST(userDataSerialization, testFile) {

    UserData u(1, 2, 3);
    {
        std::ofstream ofs("sstore.dat");
        boost::archive::text_oarchive oar(ofs);
        oar &u;
    }
    {
        std::ifstream ifs("sstore.dat");
        boost::archive::text_iarchive iar(ifs);
        UserData v;
        iar &v;
        checkUserDataEqual(u, v);
    }
}

TEST(userDataSerialization, testString) {

    UserData u(1, 2, 3);

    std::string serialized = serializeToString(u);

    checkUserDataEqual(u, deserializeFromString<UserData>(serialized));
}

TEST(userDataSerialization, testChar) {

    UserData u(1, 2, 3);
    u.p = new double[2];
    u.p[0] = 1;
    u.p[1] = 2;

    int length;
    char *buf = serializeToChar(&u, &length);

    UserData v = deserializeFromChar<UserData>(buf, length);

    delete[] buf;
    checkUserDataEqual(u, v);
}

TEST_GROUP(returnDataSerialization){void setup(){}

                          void teardown(){}};

TEST(returnDataSerialization, testString) {
    UserData u(1, 2, 3);
    Model m(1, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, AMICI_O2MODE_NONE);

    ReturnData r(&u, &m);

    std::string serialized = serializeToString(r);

    checkReturnDataEqual(r, deserializeFromString<ReturnData>(serialized));
}

