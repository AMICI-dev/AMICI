function [model] = model_nested_events_syms() 

%% CVODES OPTIONS
% set the parametrisation of the problem options are 'log', 'log10' and 'lin' (default)
model.param = 'log10';

%% STATE VARIABLES
% create state variables syms
syms Virus

% create state variable vector
x = [Virus];

%% PARAMETERS
% create parameters syms
syms V_0
syms V_0_inject
syms t_0
syms rho_V
syms delta_V

% create parameter vector
p = [V_0
     V_0_inject
     t_0
     rho_V
     delta_V];

%% CONSTANTS
% create constant vector
k = [];

%% RIGHT-HAND SIDE OF DIFFERENTIAL EQUATION
% create symbolic variable for time
syms t

% symbolic expression for right-hand side of differential equation
xdot = sym(zeros(1,1));
xdot(1) = V_0_inject*dirac(t-t_0) + Virus*rho_V*heaviside(Virus-1) - Virus*delta_V;

%% INITIAL CONDITION
x0 = sym(zeros(1,1));
x0(1) = V_0;

%% OBSERVABLE
y = sym(zeros(1,1));
y(1) = Virus;

%% EVENTS
% this part is optional and can be ommited

% events fire when there is a zero crossing of the root function
event = sym(zeros(0,1));

%% SYSTEM STRUCT
model.sym.x = x;
model.sym.k = k;
model.sym.event = event;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;

end