function [model] = model_events_syms()

% set the parametrisation of the problem options are 'log', 'log10' and
% 'lin' (default).
model.param = 'log10';
model.debug = true;

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
syms p1 p2 p3 p4

% create parameter vector 
model.sym.p = [p1,p2,p3,p4];

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
model.sym.xdot(1) = -p1*heaviside(t-p4)*x1;
% inhomogeneous
model.sym.xdot(2) = +p2*x1*exp(-0.1*t)-p3*x2 ;
model.sym.xdot(3) = -1*x3+heaviside(t-4);

%%
% INITIAL CONDITIONS

model.sym.x0 = sym(zeros(size(model.sym.x)));

model.sym.x0(1) = k1;
model.sym.x0(2) = k2;
model.sym.x0(3) = k3;

%%
% OBSERVALES

model.sym.y = sym(zeros(1,1));

model.sym.y(1) = p4 * (x1+x2+x3);


%%
% EVENTS
% this part is optional and can be ommited
syms t

% events fire when there is a -zero crossing of the root function
model.event(1) = amievent(am_ge(x2,x3),0,t);
model.event(2) = amievent(am_ge(x1,x3),0,t);

end