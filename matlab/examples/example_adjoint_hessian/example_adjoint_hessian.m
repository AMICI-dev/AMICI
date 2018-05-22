function success = example_adjoint_hessian()

%%
% COMPILATION

[exdir,~,~]=fileparts(which('example_adjoint_hessian.m'));
% compile the model
amiwrap('model_adjoint_hessian','model_adjoint_hessian_syms',exdir,1)

%%
% SIMULATION

% time vector
t = [linspace(0,4,5)];
p = [1.1,0.3,1];
k = [];

D.Y = [     1.0171
    1.3423
    1.6585
    0.9814
    0.3288];

D.Sigma_Y = 0.1*ones(size(D.Y));


options.sensi = 1;
options.sensi_meth = 'adjoint';
options.maxsteps = 1e4;
options.rtol = 1e-12;
options.atol = 1e-12;
% load mex into memory
[~] = which('simulate_model_adjoint'); % fix for inaccessability problems
sol = simulate_model_adjoint_hessian(t,log10(p),k,D,options);
g = sol.sllh;

% fprintf('First calculated gradient: \n')
% disp(g);

eps1 = [1e-7, 0, 0];
eps2 = [0, 1e-7, 0];
eps3 = [0, 0, 1e-7];

sol = simulate_model_adjoint_hessian(t,log10(p)+eps1,k,D,options);
g_p1 = sol.sllh;
sol = simulate_model_adjoint_hessian(t,log10(p)-eps1,k,D,options);
g_m1 = sol.sllh;
h1 = (g_p1-g_m1) / (2e-7);

sol = simulate_model_adjoint_hessian(t,log10(p)+eps2,k,D,options);
g_p2 = sol.sllh;
sol = simulate_model_adjoint_hessian(t,log10(p)-eps2,k,D,options);
g_m2 = sol.sllh;
h2 = (g_p2-g_m2) / (2e-7);

sol = simulate_model_adjoint_hessian(t,log10(p)+eps3,k,D,options);
g_p3 = sol.sllh;
sol = simulate_model_adjoint_hessian(t,log10(p)-eps3,k,D,options);
g_m3 = sol.sllh;
h3 = (g_p3-g_m3) / (2e-7);
H = [h1, h2, h3];

% fprintf('  Finite difference computations: \n\n');
% fprintf('  Gradient: \n');
% disp(g);
% fprintf('  Hessian: \n');
% disp(H);
% fprintf('  Hessian vector product: \n');
% disp(H * g);

fprintf('\n\n    Second order adjoint computations: \n\n')
options.sensi = 2;
options.sensi_meth = 'adjoint';
options.maxsteps = 1e4;
options.rtol = 1e-12;
options.atol = 1e-12;
% load mex into memory
[~] = which('simulate_model_adjoint'); % fix for inaccessability problems
sol1 = simulate_model_adjoint_hessian(t,log10(p),k,D,options);
% fprintf('  Likelihood: \n');
% disp(sol1.llh);
% fprintf('  Gradient: \n');
% disp(sol1.sllh);
% fprintf('  Hessian: \n');
% disp(sol1.s2llh);

if (sum(sum(abs(H - sol1.s2llh))) < 0.01)
    disp(sum(sum(abs(H - sol1.s2llh))));
    success=1;
else
    success=0;
end

end