#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.h>
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>
#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupDirac)
{
    void setup() {

    }

    void teardown() {

    }
};


amici::UserData *getTestUserData() {
    amici::UserData *udata = new amici::UserData(2,1,3);
    double par[2] = {1,3};
    double konst[1] = {1};
    udata->setParameters(par);
    udata->setConstants(konst);
    return udata;
}


/*
 * Test for mem leaks in UserData initalization / destruction
 */
TEST(groupDirac, testCreateAndFreeUserData) {
    amici::UserData *udata = getTestUserData();

    delete udata;
}

/*
 * Test for mem leaks in ExpData initalization / destruction
 */

TEST(groupDirac, testCreateAndFreeExpData) {
    amici::UserData *udata = getTestUserData();
    
    amici:: Model *model = getModel(udata);

    amici::ExpData *edata = getTestExpData(udata, model);
    
    delete model;
    delete edata;
    delete udata;
}

/*
 * Test for mem leaks in ReturnData initalization / destruction
 */

TEST(groupDirac, testCreateAndFreeReturnData) {
    amici::UserData *udata = getTestUserData();
    amici::Model *model = getModel(udata);
    amici::ReturnData *rdata = new amici::ReturnData(udata, model);
    delete model;
    delete udata;
    delete rdata;
}


TEST(groupDirac, testSimulation) {
    amici::simulateAndVerifyFromFile("/model_dirac/nosensi/");
}

TEST(groupDirac, testSimulationExpData) {

}

TEST(groupDirac, testSensitivityForward) {
   amici::simulateAndVerifyFromFile("/model_dirac/sensiforward/");
}

TEST(groupDirac, testSensitivityState) {

}

TEST(groupDirac, testSensitivityAdjoint) {

}


