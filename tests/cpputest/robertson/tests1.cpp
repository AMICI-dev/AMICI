#include "CppUTest/TestHarness.h"
#include "CppUTestExt/MockSupport.h"

#include <testfunctions.h>
#include <include/amici_hdf5.h>
#include <include/amici_interface_cpp.h>
#include <cstring>
#include "wrapfunctions.h"

TEST_GROUP(groupRobertson)
{
    void setup() {

    }

    void teardown() {

    }
};


UserData *getTestUserData() {
    Model *model = getModel();
    UserData *udata = new UserData();

    udata->qpositivex = new double[model->nx];
    udata->p = new double[model->np];
    udata->k = new double[model->nk];

    delete model;
    return udata;
}


/*
 * Test for mem leaks in UserData initalization / destruction
 */
TEST(groupRobertson, testCreateAndFreeUserData) {
    UserData *udata = getTestUserData();

    delete udata;
}

/*
 * Test for mem leaks in ExpData initalization / destruction
 */

TEST(groupRobertson, testCreateAndFreeExpData) {
    UserData udata;
    
    Model *model = getModel();

    ExpData *edata = getTestExpData(&udata, model);
    
    delete model;
    delete edata;
}

/*
 * Test for mem leaks in ReturnData initalization / destruction
 */

TEST(groupRobertson, testCreateAndFreeReturnData) {
    Model *model = getModel();

    UserData udata;
    ReturnData *rdata = new ReturnData(&udata, model);
    delete model;

    delete rdata;
}


TEST(groupRobertson, testSimulation) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_robertson/nosensi/");
    delete model;
}

TEST(groupRobertson, testSimulationExpData) {

}

TEST(groupRobertson, testSensitivityForward) {
    Model *model = getModel();
    simulateAndVerifyFromFile(model, "/model_robertson/sensiforward/");
    delete model;
}

TEST(groupRobertson, testSensitivityState) {

}

TEST(groupRobertson, testSensitivityAdjoint) {

}


