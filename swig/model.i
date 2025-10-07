%module model

// Add necessary symbols to generated header
%{
#include "amici/model.h"
using namespace amici;
%}

// remove functions that use AmiVector(Array) since that class anyways cannot
// be exposed in swig
%ignore add_adjoint_quadrature_eventUpdate;
%ignore add_adjoint_state_event_update;
%ignore add_event_objective;
%ignore add_event_objective_regularization;
%ignore add_event_objective_sensitivity;
%ignore add_observable_objective;
%ignore add_observable_objective_sensitivity;
%ignore add_partial_event_objective_sensitivity;
%ignore add_partial_observable_objective_sensitivity;
%ignore add_state_event_update;
%ignore add_state_sensitivity_event_update;
%ignore fsx_rdata;
%ignore fx_rdata;
%ignore get_adjoint_state_event_update;
%ignore get_event_time_sensitivity;
%ignore get_adjoint_state_observable_update;
%ignore get_event;
%ignore get_events;
%ignore get_event_regularization;
%ignore get_event_regularization_sensitivity;
%ignore get_event_sensitivity;
%ignore get_event_time_sensitivity;
%ignore get_explicit_roots;
%ignore get_observable;
%ignore get_observable_sensitivity;
%ignore get_expression;
%ignore init_events;
%ignore reinit_events;
%ignore initialize;
%ignore initialize_b;
%ignore initialize_state_sensitivities;
%ignore initialize_state;
%ignore reinitialize;
%ignore ModelState;
%ignore get_model_state;
%ignore set_model_state;
%ignore fx0;
%ignore fx0_fixedParameters;
%ignore fsx0;
%ignore fsx0_fixedParameters;
%ignore get_dxdotdp;
%ignore get_dxdotdp_full;
%ignore check_inite;
%ignore fJrz;
%ignore fJy;
%ignore fJz;
%ignore fdJrzdsigma;
%ignore fdJrzdz;
%ignore fdJzdsigma;
%ignore fdJzdz;
%ignore fdJydsigma;
%ignore fdeltaqB;
%ignore fdeltasx;
%ignore fdeltax;
%ignore fdeltaxB;
%ignore fdrzdp;
%ignore fdrzdx;
%ignore fdsigmaydp;
%ignore fdsigmazdp;
%ignore fdydp;
%ignore fdydx;
%ignore fdzdp;
%ignore fdzdx;
%ignore frz;
%ignore fsdx0;
%ignore fsigmay;
%ignore fsigmaz;
%ignore fsrz;
%ignore fstau;
%ignore fsz;
%ignore fw;
%ignore fy;
%ignore fz;
%ignore update_heaviside;
%ignore update_heaviside_b;
%ignore get_event_sigma;
%ignore get_event_sigma_sensitivity;
%ignore get_observable_sigma;
%ignore get_observable_sigma_sensitivity;
%ignore get_unobserved_event_sensitivity;
%ignore fdsigmaydy;
%ignore fdspline_slopesdp;
%ignore fdspline_valuesdp;
%ignore fdtotal_cldp;
%ignore fdtotal_cldx_rdata;
%ignore fdx_rdatadp;
%ignore fdx_rdatadtcl;
%ignore fdx_rdatadx_solver;
%ignore fdsigmaydy;
%ignore get_steadystate_mask_av;
%ignore initialize_splines;
%ignore initialize_spline_sensitivities;
%ignore initialize_events;
%ignore reinit_explicit_roots;
%ignore add_adjoint_quadrature_event_update;
%ignore check_finite;

%newobject amici::Model::clone;

%rename(create_solver) amici::Model::get_solver;

%extend amici::Model {
%pythoncode %{
def __deepcopy__(self, memo):
    return self.clone()
%}
};

%extend std::unique_ptr<amici::Model> {
%pythoncode %{
def __deepcopy__(self, memo):
    return self.clone()
%}
};

// Process symbols in header
%include "amici/model.h"
