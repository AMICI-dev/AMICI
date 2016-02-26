clear
[exdir,~,~]=fileparts(which('example_model_7.m'));
amiwrap('model_example_7','example_model_7_syms',exdir)

t = linspace(0,10,100);
theta = rand(8,1);

options.sensi = 1; 
options.maxsteps = 1e5; 
sol = simulate_model_example_7(t,theta,[],[],options);