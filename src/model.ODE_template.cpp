#include "TPL_MODELNAME.h"
#include <array>

namespace amici {

namespace model_TPL_MODELNAME {

std::array<const char*, TPL_NP> parameterNames = {
    TPL_PARAMETER_NAMES_INITIALIZER_LIST
};

std::array<const char*, TPL_NK> fixedParameterNames = {
    TPL_FIXED_PARAMETER_NAMES_INITIALIZER_LIST
};

std::array<const char*, TPL_NX_RDATA> stateNames = {
    TPL_STATE_NAMES_INITIALIZER_LIST
};

std::array<const char*, TPL_NY> observableNames = {
    TPL_OBSERVABLE_NAMES_INITIALIZER_LIST
};

std::array<const char*, TPL_NW> expressionNames = {
    TPL_EXPRESSION_NAMES_INITIALIZER_LIST
};

std::array<const char*, TPL_NP> parameterIds = {
    TPL_PARAMETER_IDS_INITIALIZER_LIST
};

std::array<const char*, TPL_NK> fixedParameterIds = {
    TPL_FIXED_PARAMETER_IDS_INITIALIZER_LIST
};

std::array<const char*, TPL_NX_RDATA> stateIds = {
    TPL_STATE_IDS_INITIALIZER_LIST
};

std::array<const char*, TPL_NY> observableIds = {
    TPL_OBSERVABLE_IDS_INITIALIZER_LIST
};

std::array<const char*, TPL_NW> expressionIds = {
    TPL_EXPRESSION_IDS_INITIALIZER_LIST
};

} // namespace model_TPL_MODELNAME

} // namespace amici
