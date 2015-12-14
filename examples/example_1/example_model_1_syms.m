function [model] = example_model_1_syms()

%% CVODES OPTIONS

% set the default absolute tolerance
model.atol = 1e-8; 
% set the default relative tolerance
model.rtol = 1e-8; 
% set the default maximum number of integration steps
model.maxsteps = 1e4; 
% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'log10'; 
model.debug = true;
model.recompile = false;

%% STATES

% create state syms
syms x1 x2 x3

% create state vector
x = [
x1 x2 x3
];

%% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms p1 p2 p3 p4

% create parameter vector 
p = [p1,p2,p3,p4];

%% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
syms k1 k2 k3 k4

% create parameter vector 
k = [k1 k2 k3 k4];

%% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

xdot = sym(zeros(size(x)));

% piecewise defined function
xdot(1) = -p1*heaviside(t-p4)*x1;
% inhomogeneous
xdot(2) = +p2*x1*exp(-0.1*t)-p3*x2 ;
xdot(3) = -1.5*x3;

%% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = k1;
x0(2) = k2;
x0(3) = k3;

%% OBSERVALES

y = sym(zeros(1,1));

y(1) = p4 * (x1+x2+x3);


%% EVENTS
% this part is optional and can be ommited

% events fire when there is a zero crossing of the root function
root = sym(zeros(2,1));
% x3 == x2
root(1) = x3 - x2;
% x3 == x1
root(2) = x3 - x1;

%% SYSTEM STRUCT

model.sym.x = x;
model.sym.k = k;
model.sym.root = root;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end