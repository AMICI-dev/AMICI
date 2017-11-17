#include <include/amici_model.h>
#include "wrapfunctions.h"

amici::Model *getModel(const amici::UserData *udata) {
    return new Model_model_nested_events(udata);
}

void getModelDims(int *nx, int *nk, int *np) {
    *nx = 1;
    *nk = 0;
    *np = 5;
}

