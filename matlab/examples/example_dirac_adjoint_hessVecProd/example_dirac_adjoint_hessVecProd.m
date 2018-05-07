
% COMPILATION
[exdir,~,~]=fileparts(which('example_dirac_adjoint_hessVecProd.m'));
% compile the model
amiwrap('model_dirac_adjoint_hessVecProd','model_dirac_adjoint_hessVecProd_syms',exdir,2)

%%
% SIMULATION

% time vector
tout = linspace(0,4,9);
tfine = linspace(0,4,10001);
p = [1;0.4;2;3];
k = [];
options.rtol = 1e-13;
options.atol = 1e-13;
options.maxsteps = 1e5;
options.sensi_meth = 'adjoint';

D.Y = [  0.00714742903826096
    -0.00204966058299775
    0.382159034587845
    0.33298932672138
    0.226111476113441
    0.147028440865854
    0.0882468698791813
    0.0375887796628869
    0.0373422340295005];
D.Sigma_Y = 0.01*ones(size(D.Y));


% Calculation of gradient for direction
% options.sensi = 1;
% preSol = simulate_model_dirac_adjoint(tout,log10(p),k,D,options);
% v = preSol.sllh;
% vNorm = v / sqrt(sum(v.^2));
v = [-489.2922358301243; 0; -145.7917085979483; -499.2157643879079];
vNorm = v / sqrt(sum(v.^2));

% Calculation of Hessian vector product
options.sensi = 2;
sol = simulate_model_dirac_adjoint_hessVecProd(tout,log10(p),k,D,options,v);
ASA_HVP = sol.s2llh;
format long;
fprintf('2nd order Adjoints, HVP: \n');
disp(ASA_HVP);

% Verification using Finite Differences
delta = 1e-7;
options.sensi = 1;
sol_p = simulate_model_dirac_adjoint_hessVecProd(tout, log10(p) + delta * v, k, D, options);
sol_m = simulate_model_dirac_adjoint_hessVecProd(tout, log10(p) - delta * v, k, D, options);
FD_HVP = (sol_p.sllh - sol_m.sllh)/(2*delta);

% Output of results

fprintf('Finite differences, HVP: \n');
disp(FD_HVP);

