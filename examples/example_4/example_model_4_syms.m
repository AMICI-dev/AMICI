function [model] = example_model_4_syms()

%%
% CVODES OPTIONS

model.atol = 1e-12;
model.rtol = 1e-8;
model.maxsteps = 1e4;
model.param = 'log10';
model.debug = true;

%%
% STATES

syms STAT pSTAT pSTAT_pSTAT npSTAT_npSTAT nSTAT1 nSTAT2 nSTAT3 nSTAT4 nSTAT5

x = [
STAT, pSTAT, pSTAT_pSTAT, npSTAT_npSTAT, nSTAT1, nSTAT2, nSTAT3, nSTAT4, nSTAT5 ...
];

%%
% PARAMETERS

syms p1 p2 p3 p4 init_STAT Omega_cyt Omega_nuc sp1 sp2 sp3 sp4 sp5 offset_tSTAT offset_pSTAT scale_tSTAT scale_pSTAT sigma_pSTAT sigma_tSTAT sigma_pEpoR

p = [p1,p2,p3,p4,init_STAT,sp1,sp2,sp3,sp4,sp5,offset_tSTAT,offset_pSTAT,scale_tSTAT,scale_pSTAT,sigma_pSTAT,sigma_tSTAT,sigma_pEpoR];

k = [Omega_cyt,Omega_nuc];

%%
% INPUT 
syms t
u(1) = spline_pos5(t, 0.0, sp1, 5.0, sp2, 10.0, sp3, 20.0, sp4, 60.0, sp5, 0, 0.0);

%%
% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = (Omega_nuc*p4*nSTAT5 - Omega_cyt*STAT*p1*u(1))/Omega_cyt;
xdot(2) = STAT*p1*u(1) - 2*p2*pSTAT^2;
xdot(3) = p2*pSTAT^2 - p3*pSTAT_pSTAT;
xdot(4) = -(Omega_nuc*p4*npSTAT_npSTAT - Omega_cyt*p3*pSTAT_pSTAT)/Omega_nuc;
xdot(5) = -p4*(nSTAT1 - 2*npSTAT_npSTAT);
xdot(6) = p4*(nSTAT1 - nSTAT2);
xdot(7) = p4*(nSTAT2 - nSTAT3);
xdot(8) = p4*(nSTAT3 - nSTAT4);
xdot(9) = p4*(nSTAT4 - nSTAT5);

%%
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = init_STAT;

%%
% OBSERVABLES

y = sym(zeros(3,1));

y(1) = offset_pSTAT + scale_pSTAT/init_STAT*(pSTAT + 2*pSTAT_pSTAT);
y(2) = offset_tSTAT + scale_tSTAT/init_STAT*(STAT + pSTAT + 2*(pSTAT_pSTAT));
y(3) = u(1);

%%
% SIGMA

sigma_y = sym(size(y));

sigma_y(1) = sigma_pSTAT;
sigma_y(2) = sigma_tSTAT;
sigma_y(3) = sigma_pEpoR;

%%
% SYSTEM STRUCT

model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
model.sym.sigma_y = sigma_y;
end