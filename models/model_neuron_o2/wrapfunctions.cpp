#include "amici/model.h"
#include "wrapfunctions.h"

namespace amici {

namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_model_neuron_o2::Model_model_neuron_o2());
}

} // namespace generic_model

} // namespace amici 

