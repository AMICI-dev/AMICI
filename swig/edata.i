%module edata

// Add necessary symbols to generated header
%{
#include "amici/edata.h"
using namespace amici;
%}

%ignore ConditionContext;
%ignore get_measurements_ptr;
%ignore get_noise_scales_ptr;
%ignore check_sigma_positivity;

%feature("pythonprepend") amici::ExpData::ExpData %{
    """
    Convenience wrapper for :py:class:`amici.amici.ExpData` constructors

    :param args: arguments

    :returns: ExpData Instance
    """
    if args:
        from amici.sim.sundials._numpy import ReturnDataView

        # Get the raw pointer if necessary
        if isinstance(args[0], (ExpData, ExpDataPtr, Model, ModelPtr)):
            args = (_get_ptr(args[0]), *args[1:])
        elif isinstance(args[0], ReturnDataView):
            args = (_get_ptr(args[0]["ptr"]), *args[1:])
%}

// ExpData.__repr__
%pythoncode %{
def _edata_repr(self: "ExpData"):
    n_data_y = sum(
        self.is_set_measurement(it, iy)
        for it in range(self.nt()) for
        iy in range(self.nytrue())
    )
    n_sigma_y = sum(
        self.is_set_noise_scale(it, iy)
        for it in range(self.nt())
        for iy in range(self.nytrue())
    )
    n_data_z = sum(
        self.is_set_event_measurement(ie, iz)
        for ie in range(self.nmaxevent())
        for iz in range(self.nztrue())
    )
    n_sigma_z = sum(
        self.is_set_event_noise_scale(ie, iz)
        for ie in range(self.nmaxevent())
        for iz in range(self.nztrue())
    )

    custom_simulation_settings = []
    if self.pscale:
        custom_simulation_settings.append(f"parameter scales")
    if self.fixed_parameters:
        custom_simulation_settings.append(f"fixed parameters")
    if self.fixed_parameters_pre_equilibration:
        custom_simulation_settings.append(f"pre-equilibration condition")
    if self.t_presim:
        tmp = f"pre-simulation condition (t={self.t_presim})"
        if self.fixed_parameters_presimulation:
            tmp += " with custom fixed parameters"
        custom_simulation_settings.append(tmp)
    if self.reinitialize_fixed_parameter_initial_states and self.reinitialization_state_idxs_sim:
        custom_simulation_settings.append(f"{len(self.reinitialization_state_idxs_sim)} reinitialized states (simulation)")
    if self.reinitialize_fixed_parameter_initial_states and self.reinitialization_state_idxs_presim:
        custom_simulation_settings.append(f"{len(self.reinitialization_state_idxs_presim)} reinitialized states (presimulation)")
    if self.free_parameters:
        custom_simulation_settings.append(f"free parameters")
    if self.x0:
        custom_simulation_settings.append(f"initial states")
    if self.sx0:
        custom_simulation_settings.append(f"initial state sensitivities")

    if custom_simulation_settings:
        custom_simulation_settings = " with custom " + ", ".join(custom_simulation_settings)
    else:
        custom_simulation_settings = " without custom settings"

    return "\n".join([
        self.this.__repr__()[:-1],
        f"  condition '{self.id}' starting at t={self.t_start}" + custom_simulation_settings,
        f"  {self.nt()}x{self.nytrue()} time-resolved datapoints",
        f"    ({n_data_y}/{self.nt()*self.nytrue()} measurements & {n_sigma_y}/{self.nt()*self.nytrue()} sigmas set)",
        f"  {self.nmaxevent()}x{self.nztrue()} event-resolved datapoints",
        f"    ({n_data_z}/{self.nmaxevent()*self.nztrue()} measurements & {n_sigma_z}/{self.nmaxevent()*self.nztrue()} sigmas set)",
        ">"
    ])
%}
%extend amici::ExpData {
%pythoncode %{
def __repr__(self):
    return _edata_repr(self)

def __eq__(self, other):
    return other.__class__ == self.__class__ and __eq__(self, other)

def __deepcopy__(self, memo):
    # invoke copy constructor
    return type(self)(self)

def __reduce__(self):
    from amici.sim.sundials._swig_wrappers import restore_edata

    return (
        restore_edata,
        (
            # ExpData ctor arguments
            (
                self.nytrue(),
                self.nztrue(),
                self.nmaxevent(),
                self.get_timepoints(),
                self.get_measurements(),
                self.get_noise_scales(),
                self.get_event_measurements(),
                self.get_event_noise_scales(),
            ),
            dict(self)
        ),
        {}
    )
%}
};
%extend std::unique_ptr<amici::ExpData> {
%pythoncode %{
def __repr__(self):
    return _edata_repr(self)
%}
};


// Process symbols in header
%include "amici/edata.h"
