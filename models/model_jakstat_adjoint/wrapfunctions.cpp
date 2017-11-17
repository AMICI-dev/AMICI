#include <include/amici_model.h>
#include "wrapfunctions.h"

amici::Model *getModel(const amici::UserData *udata) {
    return new Model_model_jakstat_adjoint(udata);
}

void getModelDims(int *nx, int *nk, int *np) {
    *nx = 9;
    *nk = 2;
    *np = 17;
}

