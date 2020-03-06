%module model

// Add necessary symbols to generated header
%{
#include "amici/model.h"
using namespace amici;
%}

// remove functions that use AmiVector(Array) since that class anyways cannot
// be exposed in swig
%ignore addAdjointQuadratureEventUpdate;
%ignore addAdjointStateEventUpdate;
%ignore addEventObjective;
%ignore addEventObjectiveRegularization;
%ignore addEventObjectiveSensitivity;
%ignore addObservableObjective;
%ignore addObservableObjectiveSensitivity;
%ignore addPartialEventObjectiveSensitivity;
%ignore addPartialObservableObjectiveSensitivity;
%ignore addStateEventUpdate;
%ignore addStateSensitivityEventUpdate;
%ignore fsx_rdata;
%ignore fx_rdata;
%ignore getAdjointStateEventUpdate;
%ignore getEventTimeSensitivity;
%ignore getAdjointStateObservableUpdate;
%ignore getEvent;
%ignore getEventRegularization;
%ignore getEventRegularizationSensitivity;
%ignore getEventSensitivity;
%ignore getEventTimeSensitivity;
%ignore getObservable;
%ignore getObservableSensitivity;
%ignore getExpression;
%ignore initHeaviside;
%ignore initialize;
%ignore initializeB;
%ignore initializeStateSensitivities;
%ignore initializeStates;
%ignore dxdotdp;
%ignore dxdotdp_implicit;
%ignore dxdotdp_explicit;

// Process symbols in header
%include "amici/model.h"
