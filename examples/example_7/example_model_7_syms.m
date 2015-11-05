function [model] = example_model_7_syms()

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;
model.debug = true;

% STATES

syms p_y_1_0 p_y_0_1 mu_1_y_1_0 mu_2_y_1_0 C_1_1_y_1_0 C_1_2_y_1_0 C_2_2_y_1_0 mu_1_y_0_1 mu_2_y_0_1 C_1_1_y_0_1 C_1_2_y_0_1 C_2_2_y_0_1

x = [
p_y_1_0, p_y_0_1, mu_1_y_1_0, mu_2_y_1_0, C_1_1_y_1_0, C_1_2_y_1_0, C_2_2_y_1_0, mu_1_y_0_1, mu_2_y_0_1, C_1_1_y_0_1, C_1_2_y_0_1, C_2_2_y_0_1 ...
];

% PARAMETERS

syms tau_on tau_off k_m gamma_m k_p gamma_p tau_on_p Omega 

p = [tau_on,tau_off,k_m,gamma_m,k_p,gamma_p,tau_on_p,Omega];

% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

f = sym(zeros(size(x)));

f(1) = p_y_0_1*tau_off - p_y_1_0*tau_on - mu_2_y_1_0*p_y_1_0*tau_on_p;
f(2) = p_y_1_0*tau_on - p_y_0_1*tau_off + mu_2_y_1_0*p_y_1_0*tau_on_p;
f(3) = p_y_0_1*(mu_1_y_0_1*tau_off - mu_1_y_1_0*tau_off) - C_1_2_y_1_0*p_y_1_0*tau_on_p - gamma_m*mu_1_y_1_0*p_y_1_0;
f(4) = p_y_0_1*(mu_2_y_0_1*tau_off - mu_2_y_1_0*tau_off) - C_2_2_y_1_0*p_y_1_0*tau_on_p - gamma_p*mu_2_y_1_0*p_y_1_0 + k_p*mu_1_y_1_0*p_y_1_0;
f(5) = p_y_0_1*(C_1_1_y_0_1*tau_off - C_1_1_y_1_0*tau_off + mu_1_y_0_1^2*tau_off + mu_1_y_1_0^2*tau_off - 2*mu_1_y_0_1*mu_1_y_1_0*tau_off) - 2*C_1_1_y_1_0*gamma_m*p_y_1_0 + gamma_m*mu_1_y_1_0*p_y_1_0;
f(6) = C_1_1_y_1_0*k_p*p_y_1_0 - C_1_2_y_1_0*gamma_p*p_y_1_0 - C_1_2_y_1_0*gamma_m*p_y_1_0 + C_1_2_y_0_1*p_y_0_1*tau_off - C_1_2_y_1_0*p_y_0_1*tau_off + mu_1_y_0_1*mu_2_y_0_1*p_y_0_1*tau_off - mu_1_y_0_1*mu_2_y_1_0*p_y_0_1*tau_off - mu_1_y_1_0*mu_2_y_0_1*p_y_0_1*tau_off + mu_1_y_1_0*mu_2_y_1_0*p_y_0_1*tau_off;
f(7) = p_y_0_1*(C_2_2_y_0_1*tau_off - C_2_2_y_1_0*tau_off + mu_2_y_0_1^2*tau_off + mu_2_y_1_0^2*tau_off - 2*mu_2_y_0_1*mu_2_y_1_0*tau_off) - 2*C_2_2_y_1_0*gamma_p*p_y_1_0 + 2*C_1_2_y_1_0*k_p*p_y_1_0 + gamma_p*mu_2_y_1_0*p_y_1_0 + k_p*mu_1_y_1_0*p_y_1_0;
f(8) = p_y_0_1*(k_m - gamma_m*mu_1_y_0_1) + C_1_2_y_1_0*p_y_1_0*tau_on_p - mu_1_y_0_1*p_y_1_0*tau_on + mu_1_y_1_0*p_y_1_0*tau_on - mu_1_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p + mu_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p;
f(9) = mu_2_y_1_0^2*p_y_1_0*tau_on_p - p_y_0_1*(gamma_p*mu_2_y_0_1 - k_p*mu_1_y_0_1) + C_2_2_y_1_0*p_y_1_0*tau_on_p - mu_2_y_0_1*p_y_1_0*tau_on + mu_2_y_1_0*p_y_1_0*tau_on - mu_2_y_0_1*mu_2_y_1_0*p_y_1_0*tau_on_p;
f(10) = p_y_0_1*(k_m - 2*C_1_1_y_0_1*gamma_m + gamma_m*mu_1_y_0_1) - C_1_1_y_0_1*p_y_1_0*(tau_on + mu_2_y_1_0*tau_on_p) + p_y_1_0*tau_on*(mu_1_y_0_1 - mu_1_y_1_0)^2 + C_1_1_y_1_0*p_y_1_0*tau_on - C_1_2_y_1_0*p_y_1_0*tau_on_p*(2*mu_1_y_0_1 - 2*mu_1_y_1_0) + mu_2_y_1_0*p_y_1_0*tau_on_p*(mu_1_y_0_1 - mu_1_y_1_0)^2 + C_1_1_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p;
f(11) = C_1_2_y_1_0*p_y_1_0*tau_on - C_1_2_y_0_1*p_y_1_0*(tau_on + mu_2_y_1_0*tau_on_p) - p_y_0_1*(C_1_2_y_0_1*gamma_m + C_1_2_y_0_1*gamma_p - C_1_1_y_0_1*k_p) - C_1_2_y_1_0*p_y_1_0*tau_on_p*(mu_2_y_0_1 - mu_2_y_1_0) - C_2_2_y_1_0*p_y_1_0*tau_on_p*(mu_1_y_0_1 - mu_1_y_1_0) + p_y_1_0*tau_on*(mu_1_y_0_1 - mu_1_y_1_0)*(mu_2_y_0_1 - mu_2_y_1_0) + C_1_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p + mu_2_y_1_0*p_y_1_0*tau_on_p*(mu_1_y_0_1 - mu_1_y_1_0)*(mu_2_y_0_1 - mu_2_y_1_0);
f(12) = p_y_0_1*(2*C_1_2_y_0_1*k_p - 2*C_2_2_y_0_1*gamma_p + gamma_p*mu_2_y_0_1 + k_p*mu_1_y_0_1) - C_2_2_y_0_1*p_y_1_0*(tau_on + mu_2_y_1_0*tau_on_p) + p_y_1_0*tau_on*(mu_2_y_0_1 - mu_2_y_1_0)^2 + C_2_2_y_1_0*p_y_1_0*tau_on - C_2_2_y_1_0*p_y_1_0*tau_on_p*(2*mu_2_y_0_1 - 2*mu_2_y_1_0) + mu_2_y_1_0*p_y_1_0*tau_on_p*(mu_2_y_0_1 - mu_2_y_1_0)^2 + C_2_2_y_1_0*mu_2_y_1_0*p_y_1_0*tau_on_p;
M = diag(sym([[1], [1], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1]]));
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = 1;
x0(2) = 1/10000000000;
x0(3) = 4;
x0(4) = 10;
x0(5) = 1/10000000000;
x0(6) = 1/10000000000;
x0(7) = 1/10000000000;
x0(8) = 4;
x0(9) = 10;
x0(10) = 1/10000000000;
x0(11) = 1/10000000000;
x0(12) = 1/10000000000;
dx0 = sym(zeros(size(x)));
dx0(1) = tau_off/10000000000 - tau_on - 10*tau_on_p;
dx0(2) = tau_on - tau_off/10000000000 + 10*tau_on_p;
dx0(3) = - 4*gamma_m - tau_on_p/10000000000;
dx0(4) = 4*k_p - 10*gamma_p - tau_on_p/10000000000;
dx0(5) = (19999999999*gamma_m)/5000000000;
dx0(6) = k_p/10000000000 - gamma_p/10000000000 - gamma_m/10000000000;
dx0(7) = (49999999999*gamma_p)/5000000000 + (20000000001*k_p)/5000000000;
dx0(8) = k_m - 4*gamma_m + tau_on_p;
dx0(9) = 4*k_p - 10*gamma_p + tau_on_p;
dx0(10) = (19999999999*gamma_m)/5000000000 + k_m;
dx0(11) = k_p/10000000000 - gamma_p/10000000000 - gamma_m/10000000000;
dx0(12) = (49999999999*gamma_p)/5000000000 + (20000000001*k_p)/5000000000;

% OBSERVABLES

y = sym(zeros(14,1));

y(1) = p_y_1_0;
y(2) = p_y_0_1;
y(3) = mu_1_y_0_1*p_y_0_1 + mu_1_y_1_0*p_y_1_0;
y(4) = mu_2_y_0_1*p_y_0_1 + mu_2_y_1_0*p_y_1_0;
y(5) = p_y_0_1*p_y_1_0^2 + p_y_1_0*(p_y_1_0 - 1)^2;
y(6) = p_y_0_1*p_y_1_0*(p_y_0_1 + p_y_1_0 - 2);
y(7) = p_y_1_0*(p_y_1_0 - 1)*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0);
y(8) = p_y_1_0*(p_y_1_0 - 1)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0);
y(9) = p_y_0_1^2*p_y_1_0 + p_y_0_1*(p_y_0_1 - 1)^2;
y(10) = p_y_0_1*(p_y_0_1 - 1)*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0);
y(11) = p_y_0_1*(p_y_0_1 - 1)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0);
y(12) = p_y_0_1*(C_1_1_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_1_1_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)^2);
y(13) = p_y_0_1*(C_1_2_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)) + p_y_1_0*(C_1_2_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0));
y(14) = p_y_0_1*(C_2_2_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_2_2_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^2);

% SYSTEM STRUCT

model.sym.x = x;
model.sym.u = u;
model.sym.f = f;
model.sym.M = M;
model.sym.p = p;
model.sym.x0 = x0;
model.sym.dx0 = dx0;
model.sym.y = y;
end