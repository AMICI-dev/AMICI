function [model] = model_robertson_syms()


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
syms p1 p2 p3

% create parameter vector 
model.sym.p = [p1,p2,p3];

% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'log10'; 


%%
% CONSTANTS ( for these no sensitivities will be computed )
% this part is optional and can be ommited

% create parameter syms
syms k1

% create parameter vector 
model.sym.k = [k1];

%%
% SYSTEM EQUATIONS

% create symbolic variable for time
syms t

model.sym.xdot = sym(zeros(size(model.sym.x)));
model.sym.xdot(1) = -p1*x1+p2*x2*x3;
model.sym.xdot(2) = p1*x1-p2*x2*x3-p3*x2^2;
model.sym.xdot(3) = x1+x2+x3-1;

model.sym.M = [1 0 0; 0 1 0; 0 0 0];

%%
% INITIAL CONDITIONS

model.sym.x0 = sym(zeros(size(model.sym.x)));

model.sym.x0(1) = k1;
model.sym.x0(2) = 0;
model.sym.x0(3) = 0;

model.sym.dx0 = sym(zeros(size(model.sym.x)));

%%
% OBSERVALES

model.sym.y = model.sym.x;
model.sym.y(2) = 1e4*model.sym.x(2);
end