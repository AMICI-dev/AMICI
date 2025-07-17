#include "wrapfunctions.h"
#include "model_dirac_py.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_model_dirac_py::Model_model_dirac_py()
    );
}

} // namespace generic_model

} // namespace amici
