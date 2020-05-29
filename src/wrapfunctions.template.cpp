#include "amici/model.h"
#include "wrapfunctions.h"
#include "TPL_MODELNAME.h"

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(new Model_TPL_MODELNAME());
}
