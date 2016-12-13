function [model] = model_adjoint_syms()

%%
% STATES

% create state syms
syms x1

% create state vector
model.sym.x = [x1];

%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms p1 p2 p3 

% create parameter vector 
model.sym.p = [p1 p2 p3];


% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'log10';

%%
% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

model.sym.xdot = sym(zeros(size(model.sym.x)));

% piecewise defined function
model.sym.xdot(1) = -p1*x1*heaviside(t-2) + p2;

%%
% INITIAL CONDITIONS

model.sym.x0 = sym(zeros(size(model.sym.x)));

model.sym.x0(1) = p3;

%%
% OBSERVALES

model.sym.y = sym(zeros(1,1));

model.sym.y(1) = x1;

end