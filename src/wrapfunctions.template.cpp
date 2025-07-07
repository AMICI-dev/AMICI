#include "wrapfunctions.h"
#include "TPL_MODELNAME.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_TPL_MODELNAME::Model_TPL_MODELNAME()
    );
}

} // namespace generic_model

} // namespace amici
