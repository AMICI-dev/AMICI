% function [model] = MEC_2_LD_2_c_f_AMICI_genExp_syms(f0_user)
function [model] = geneExpression_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;
model.adjoint = false;
model.forward = false;
model.debug = true;

% STATES

syms mu_1 mu_2 mu_3 mu_4 C_1_1 C_1_2 C_1_3 C_1_4 C_2_2 C_2_3 C_2_4 C_3_3 C_3_4 C_4_4

x = [
mu_1, mu_2, mu_3, mu_4, C_1_1, C_1_2, C_1_3, C_1_4, C_2_2, C_2_3, C_2_4, C_3_3, C_3_4, C_4_4 ...
];

% PARAMETERS

syms tau_on tau_off k_m gamma_m k_p gamma_p tau_on_p scaleP offsetP r0 

% KAPPA (constant parameters)

syms Omega indmu1 indmu2 indmu3 indmu4 indC1 indC2 indC3 indC4 indC5 indC6 indC7 indC8 indC9 indC10 kmu01 kmu02 kmu03 kmu04 kC01 kC02 kC03 kC04 kC05 kC06 kC07 kC08 kC09 kC010 

syms t

p = [tau_on,tau_off,k_m,gamma_m,k_p,gamma_p,tau_on_p,scaleP,offsetP,r0];

k = [Omega,indmu1,indmu2,indmu3,indmu4,indC1,indC2,indC3,indC4,indC5,indC6,indC7,indC8,indC9,indC10,kmu01,kmu02,kmu03,kmu04,kC01,kC02,kC03,kC04,kC05,kC06,kC07,kC08,kC09,kC010];

if nargin > 0
   f0_user = varargin{1};
   if ~isnumeric(f0_user)
      p_user = setdiff(symvar(f0_user),p);
      % ADDITIONAL PARAMETERS IN INITIAL CONDITIONS
      p = [p,p_user];
   end
	fmu01 = f0_user(1); 
	fmu02 = f0_user(2); 
	fmu03 = f0_user(3); 
	fmu04 = f0_user(4); 
	fC01 = f0_user(5); 
	fC02 = f0_user(6); 
	fC03 = f0_user(7); 
	fC04 = f0_user(8); 
	fC05 = f0_user(9); 
	fC06 = f0_user(10); 
	fC07 = f0_user(11); 
	fC08 = f0_user(12); 
	fC09 = f0_user(13); 
	fC010 = f0_user(14); 
else
	fmu01 = 1; 
	fmu02 = 0; 
	fmu03 = r0; 
	fmu04 = 0; 
	fC01 = 0; 
	fC02 = 0; 
	fC03 = 0; 
	fC04 = 0; 
	fC05 = 0; 
	fC06 = 0; 
	fC07 = 0; 
	fC08 = 0; 
	fC09 = 0; 
	fC010 = 0; 
end
% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = -(Omega*mu_1*tau_on - Omega*mu_2*tau_off + C_1_4*Omega^2*tau_on_p + Omega^2*mu_1*mu_4*tau_on_p)/Omega;
xdot(2) = (Omega*mu_1*tau_on - Omega*mu_2*tau_off + C_1_4*Omega^2*tau_on_p + Omega^2*mu_1*mu_4*tau_on_p)/Omega;
xdot(3) = -(Omega*gamma_m*mu_3 - Omega*k_m*mu_2)/Omega;
xdot(4) = -(Omega*gamma_p*mu_4 - Omega*k_p*mu_3)/Omega;
xdot(5) = (Omega*mu_1*tau_on + Omega*mu_2*tau_off - 2*C_1_1*Omega^2*tau_on + 2*C_1_2*Omega^2*tau_off + C_1_4*Omega^2*tau_on_p - 2*C_1_1*Omega^3*mu_4*tau_on_p - 2*C_1_4*Omega^3*mu_1*tau_on_p + Omega^2*mu_1*mu_4*tau_on_p)/Omega^2;
xdot(6) = -(Omega*mu_1*tau_on + Omega*mu_2*tau_off - C_1_1*Omega^2*tau_on + C_1_2*Omega^2*tau_on + C_1_2*Omega^2*tau_off - C_2_2*Omega^2*tau_off + C_1_4*Omega^2*tau_on_p - C_1_1*Omega^3*mu_4*tau_on_p - C_1_4*Omega^3*mu_1*tau_on_p + C_1_2*Omega^3*mu_4*tau_on_p + C_2_4*Omega^3*mu_1*tau_on_p + Omega^2*mu_1*mu_4*tau_on_p)/Omega^2;
xdot(7) = -(C_1_3*Omega^2*gamma_m - C_1_2*Omega^2*k_m + C_1_3*Omega^2*tau_on - C_2_3*Omega^2*tau_off + C_1_3*Omega^3*mu_4*tau_on_p + C_3_4*Omega^3*mu_1*tau_on_p)/Omega^2;
xdot(8) = -(C_1_4*Omega^2*gamma_p - C_1_3*Omega^2*k_p + C_1_4*Omega^2*tau_on - C_2_4*Omega^2*tau_off + C_1_4*Omega^3*mu_4*tau_on_p + C_4_4*Omega^3*mu_1*tau_on_p)/Omega^2;
xdot(9) = (Omega*mu_1*tau_on + Omega*mu_2*tau_off + 2*C_1_2*Omega^2*tau_on - 2*C_2_2*Omega^2*tau_off + C_1_4*Omega^2*tau_on_p + 2*C_1_2*Omega^3*mu_4*tau_on_p + 2*C_2_4*Omega^3*mu_1*tau_on_p + Omega^2*mu_1*mu_4*tau_on_p)/Omega^2;
xdot(10) = (C_2_2*Omega^2*k_m - C_2_3*Omega^2*gamma_m + C_1_3*Omega^2*tau_on - C_2_3*Omega^2*tau_off + C_1_3*Omega^3*mu_4*tau_on_p + C_3_4*Omega^3*mu_1*tau_on_p)/Omega^2;
xdot(11) = (C_2_3*Omega^2*k_p - C_2_4*Omega^2*gamma_p + C_1_4*Omega^2*tau_on - C_2_4*Omega^2*tau_off + C_1_4*Omega^3*mu_4*tau_on_p + C_4_4*Omega^3*mu_1*tau_on_p)/Omega^2;
xdot(12) = (Omega*gamma_m*mu_3 + Omega*k_m*mu_2 - 2*C_3_3*Omega^2*gamma_m + 2*C_2_3*Omega^2*k_m)/Omega^2;
xdot(13) = -(C_3_4*Omega^2*gamma_m + C_3_4*Omega^2*gamma_p - C_2_4*Omega^2*k_m - C_3_3*Omega^2*k_p)/Omega^2;
xdot(14) = (Omega*gamma_p*mu_4 + Omega*k_p*mu_3 - 2*C_4_4*Omega^2*gamma_p + 2*C_3_4*Omega^2*k_p)/Omega^2;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = (indmu1*kmu01 - fmu01*(indmu1 - 1))/Omega;
x0(2) = (indmu2*kmu02 - fmu02*(indmu2 - 1))/Omega;
x0(3) = (indmu3*kmu03 - fmu03*(indmu3 - 1))/Omega;
x0(4) = (indmu4*kmu04 - fmu04*(indmu4 - 1))/Omega;
x0(5) = (indC1*kC01 - fC01*(indC1 - 1))/Omega^2;
x0(6) = (indC2*kC02 - fC02*(indC2 - 1))/Omega^2;
x0(7) = (indC3*kC03 - fC03*(indC3 - 1))/Omega^2;
x0(8) = (indC4*kC04 - fC04*(indC4 - 1))/Omega^2;
x0(9) = (indC5*kC05 - fC05*(indC5 - 1))/Omega^2;
x0(10) = (indC6*kC06 - fC06*(indC6 - 1))/Omega^2;
x0(11) = (indC7*kC07 - fC07*(indC7 - 1))/Omega^2;
x0(12) = (indC8*kC08 - fC08*(indC8 - 1))/Omega^2;
x0(13) = (indC9*kC09 - fC09*(indC9 - 1))/Omega^2;
x0(14) = (indC10*kC010 - fC010*(indC10 - 1))/Omega^2;

% OBSERVABLES

y = sym(zeros(2,1));

y(1) = offsetP + mu_4*scaleP;
y(2) = C_4_4*scaleP^2;

% SYSTEM STRUCT

model.sym.nmx = 0;
model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
% Additional fields for the prespecified length of kappa
model.sym.nk1 = 1;
end