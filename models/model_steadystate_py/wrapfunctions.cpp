#include "wrapfunctions.h"
#include "model_steadystate_py.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> get_model() {
    return std::unique_ptr<amici::Model>(
        new amici::model_model_steadystate_py::Model_model_steadystate_py()
    );
}

} // namespace generic_model

} // namespace amici
