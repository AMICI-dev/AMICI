#include "wrapfunctions.h"
#include "model_jakstat_adjoint_py.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> get_model() {
    return std::unique_ptr<amici::Model>(
        new amici::model_model_jakstat_adjoint_py::Model_model_jakstat_adjoint_py()
    );
}

} // namespace generic_model

} // namespace amici
