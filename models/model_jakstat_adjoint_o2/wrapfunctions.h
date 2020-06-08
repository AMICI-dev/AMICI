#ifndef _amici_wrapfunctions_h
#define _amici_wrapfunctions_h

#include "model_jakstat_adjoint_o2.h"

namespace amici {

namespace generic_model {

std::unique_ptr<amici::Model> getModel();

} // namespace generic_model

} // namespace amici 

#endif /* _amici_wrapfunctions_h */
