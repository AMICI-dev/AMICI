function [model] = example_model_5_syms()

%%
% CVODES OPTIONS

% set the default absolute tolerance
model.atol = 1e-8; 
% set the default relative tolerance
model.rtol = 1e-8; 
% set the default maximum number of integration steps
model.maxsteps = 1e4; 
% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'log10';

%%
% STATES

% create state syms
syms x1 x2

% create state vector
x = [ x1 x2 ];

%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms p1 p2 p3 p4

% create parameter vector 
p = [p1,p2,p3,p4];

%%
% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

xdot = sym(zeros(size(x)));

% piecewise defined function
xdot(1) = -p1*x1 + dirac(t-p2);
% inhomogeneous
xdot(2) = p3*x1 - p4*x2 ;

%%
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = 0;
x0(2) = 0;

%%
% OBSERVALES

y = sym(zeros(1,1));

y(1) = x2;

%%
% SYSTEM STRUCT

model.sym.x = x;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.y = y;
end