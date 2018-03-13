#include "amici/amici_model.h"
#include "wrapfunctions.h"

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(new Model_model_nested_events());
}

