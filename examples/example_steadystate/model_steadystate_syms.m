function [model] = model_steadystate_syms()


%%
% STATES

% create state syms
syms x1 x2 x3

% create state vector
model.sym.x = [
x1 x2 x3
];

%%
% PARAMETERS ( for these sensitivities will be computed )

% create parameter syms
syms p1 p2 p3 p4 p5

% create parameter vector 
model.sym.p = [p1,p2,p3,p4,p5];

% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'log10'; 


%%
% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
syms k1 k2 k3 k4

% create parameter vector 
model.sym.k = [k1 k2 k3 k4];

%%
% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

model.sym.xdot = sym(zeros(size(model.sym.x)));

% piecewise defined function
model.sym.xdot(1) = -2*p1*x1^2 - p2*x1*x2 + 2*p3*x2 + p4*x3 + p5;
% inhomogeneous
model.sym.xdot(2) = +p1*x1^2 - p2*x1*x2 - p3*x2 + p4*x3;
model.sym.xdot(3) = p2*x1*x2 - p4*x3 - k4*x3;

%%
% INITIAL CONDITIONS

model.sym.x0 = sym(zeros(size(model.sym.x)));

model.sym.x0(1) = k1;
model.sym.x0(2) = k2;
model.sym.x0(3) = k3;

%%
% OBSERVALES

model.sym.y = model.sym.x;
end