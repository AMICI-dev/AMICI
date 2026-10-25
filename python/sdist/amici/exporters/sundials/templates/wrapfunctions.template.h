#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h

#include <memory>

#include "amici/model.h"

namespace amici::generic_model {

/**
 * @brief Wrapper function to instantiate the linked Amici model without knowing
 * the name at compile time.
 * @return Model instance
 */
std::unique_ptr<amici::Model> get_model();

} // namespace amici::generic_model


#endif /* _amici_wrapfunctions_h */
