amiwrap('geneExpression','geneExpression_syms',pwd)

t = linspace(0,100,500);
theta = [0.3;0.3;10;1;4;1;0.015;1;0;1e-10];
kappa = [1;zeros(14,1);zeros(14,1)];

sol = simulate_geneExpression(t,theta,kappa,[]);