#include <include/amici_model.h>
#include "wrapfunctions.h"

amici::Model *getModel(const amici::UserData *udata) {
    return new Model_model_dirac(udata);
}

void getModelDims(int *nx, int *nk, int *np) {
    *nx = 2;
    *nk = 0;
    *np = 4;
}

