#include <amici/defines.h>
#include <array>
#include <string_view>

namespace amici::model_TPL_MODELNAME {

// clang-format off

extern const std::array<std::string_view const, TPL_NP> free_parameter_names = {
    TPL_FREE_PARAMETER_NAMES_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NK> fixed_parameter_names = {
    TPL_FIXED_PARAMETER_NAMES_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NX_RDATA> state_names = {
    TPL_STATE_NAMES_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NX_SOLVER> state_names_solver = {
    TPL_STATE_NAMES_SOLVER_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NY> observable_names = {
    TPL_OBSERVABLE_NAMES_INITIALIZER_LIST
};

std::array<const ObservableScaling, TPL_NY> observable_scalings = {
    TPL_OBSERVABLE_TRAFO_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NW> expression_names = {
    TPL_EXPRESSION_NAMES_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NP> free_parameter_ids = {
    TPL_FREE_PARAMETER_IDS_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NK> fixed_parameter_ids = {
    TPL_FIXED_PARAMETER_IDS_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NX_RDATA> state_ids = {
    TPL_STATE_IDS_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NX_SOLVER> state_ids_solver = {
    TPL_STATE_IDS_SOLVER_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NY> observable_ids = {
    TPL_OBSERVABLE_IDS_INITIALIZER_LIST
};

extern const std::array<std::string_view const, TPL_NW> expression_ids = {
    TPL_EXPRESSION_IDS_INITIALIZER_LIST
};

std::array<int, TPL_NX_SOLVER> state_idxs_solver = {
    TPL_STATE_IDXS_SOLVER_INITIALIZER_LIST
};

// clang-format on

} // namespace amici::model_TPL_MODELNAME
