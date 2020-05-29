#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h

#include <memory>

#include "amici/model.h"

/**
 * @brief Wrapper function to instantiate the linked Amici model without knowing
 * the name at compile time.
 * @return Model instance
 */
std::unique_ptr<amici::Model> getModel();

#endif /* _amici_wrapfunctions_h */
