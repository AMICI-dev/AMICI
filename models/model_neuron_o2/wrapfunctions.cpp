#include <include/amici_model.h>
#include "wrapfunctions.h"

std::unique_ptr<amici::Model> getModel(const amici::UserData *udata) {
    return std::unique_ptr<amici::Model>(new Model_model_neuron_o2(udata));
}

void getModelDims(int *nx, int *nk, int *np) {
    if(nx)
        *nx = 10;
    if(nk)
        *nk = 2;
    if(np)
        *np = 4;
}

