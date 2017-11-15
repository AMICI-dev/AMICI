#include <include/amici_model.h>
#include "wrapfunctions.h"

amici::Model *getModel(const amici::UserData *udata) {
    return new Model_model_events(udata);
}

void getModelDims(int *nx, int *nk, int *np) {
    *nx = 3;
    *nk = 4;
    *np = 4;
}

