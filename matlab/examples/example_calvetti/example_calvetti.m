function example_calvetti()
%%
% COMPILATION

[exdir,~,~]=fileparts(which('example_calvetti.m'));
% compile the model
amiwrap('model_calvetti','model_calvetti_syms',exdir)

% time vector
t = linspace(0,20,201);
p = [];
k = [0.29 0.74 0.44 0.08 0.27 0.18];

options = amioption('sensi',0,...
    'maxsteps',1e4);
D = amidata(length(t),6,0,2,4);
% load mex into memory
[~] = which('simulate_model_calvetti'); % fix for inaccessability problems
sol = simulate_model_calvetti(t,p,k,D,options);

tic
sol = simulate_model_calvetti(t,p,k,D,options);
disp(['Time elapsed with cvodes: ' num2str(toc) ])

%%
% ODE15S

y0 = [k(1); k(3); k(5); 1; 1; 1;];
M = [1 0 0 0 0 0 
     0 1 0 0 0 0
     0 0 1 0 0 0
     0 0 0 0 0 0
     0 0 0 0 0 0
     0 0 0 0 0 0];
 
function [xdot] = dae_system(t,x,p,k,it)
    if it<3
        h0 = 0;
    else
        h0 = 1;
    end
    if it<4
        h1 = 0;
    else
        h1 = 1;
    end
    xdot = [
    h0/31 - h1/31 - (100*x(1))/(899*k(1)) + (2*x(1)*(k(2)/2 - 1))/(31*k(1)*(x(4)*((k(2)*k(1)^2)/x(1)^2 + (k(4)*k(3)^2)/x(2)^2) + x(5)*((k(4)*k(3)^2)/x(2)^2 + (k(6)*k(5)^2)/x(3)^2) + (k(6)*k(5)^2*x(6))/x(3)^2)) + 129/899;
    (2*x(2)*(k(2) + k(4)/2 - 1))/(163*k(3)*(x(5)*((k(4)*k(3)^2)/x(2)^2 + (k(6)*k(5)^2)/x(3)^2) + (k(6)*k(5)^2*x(6))/x(3)^2)) - (100*x(2))/(8313*k(3)) + 151/8313;
    (x(3)^3*(k(2) + k(4) +k(6)/2 - 1))/(61*k(6)*k(5)^3*x(6)) - x(3)/(121999878*k(5)) + 500000/60999939;
    h1/31 - h0/31 - x(4) + (100*x(1))/(899*k(1)) + (2*x(1)^2)/(k(2)*k(1)^2) - (x(1)^2*(x(4)*((k(2)*k(1)^2)/x(1)^2 + (k(4)*k(3)^2)/x(2)^2) + x(5)*((k(4)*k(3)^2)/x(2)^2 + (k(6)*k(5)^2)/x(3)^2) + (k(6)*k(5)^2*x(6))/x(3)^2))/(k(2)*k(1)^2) - (2*x(1)*(k(2)/2 - 1))/(31*k(1)*(x(4)*((k(2)*k(1)^2)/x(1)^2 + (k(4)*k(3)^2)/x(2)^2) + x(5)*((k(4)*k(3)^2)/x(2)^2 + (k(6)*k(5)^2)/x(3)^2) + (k(6)*k(5)^2*x(6))/x(3)^2)) - 129/899;
    x(4) - x(5) + (100*x(2))/(8313*k(3)) - (2*x(2)*(k(2) + k(4)/2 - 1))/(163*k(3)*(x(5)*((k(4)*k(3)^2)/x(2)^2 + (k(6)*k(5)^2)/x(3)^2) + (k(6)*k(5)^2*x(6))/x(3)^2)) - 151/8313;
    x(5) - x(6) + x(3)/(121999878*k(5)) - (x(3)^3*(k(2) + k(4) +k(6)/2 - 1))/(61*k(6)*k(5)^3*x(6)) - 500000/60999939;
];
end

options_ode15s = odeset('Mass',M, ...
                        'RelTol',options.rtol, ...
                        'AbsTol',options.atol);

tic
X_ode15s = [];
tt = t;
tstop = [0,10,12,20];
for it = 2:4
    [~, y] = ode15s(@(t,x) dae_system(t,x,p,k,it),t(t>=tstop(it-1) & t<=tstop(it)),y0,options_ode15s);
    X_ode15s = [X_ode15s;y(1:end-1,:)];
    y0 = y(end,:);
end
X_ode15s = [X_ode15s;y(end,:)];
disp(['Time elapsed with ode15s: ' num2str(toc) ])

%%
% PLOTTING
if(usejava('jvm'))
    figure
    c_x = get(gca,'ColorOrder');
    subplot(2,1,1)
    for ix = 1:size(sol.x,2)
        plot(t,sol.x(:,ix),'.-','Color',c_x(ix,:))
        hold on
        plot(t,X_ode15s(:,ix),'d','Color',c_x(ix,:))
    end
    legend('x1','x1_{ode15s}','x2','x2_{ode15s}', ...
           'x3','x3_{ode15s}','x4','x4_{ode15s}', ...
           'x5','x5_{ode15s}','x6','x6_{ode15s}', ...
           'Location','NorthEastOutside')
    legend boxoff
    xlabel('time t')
    ylabel('x')
    box on
    subplot(2,1,2)
    plot(t,abs(sol.x-X_ode15s),'--')
    set(gca,'YScale','log')
    legend('error x1','error x2','error x3','error x4','error x5','error x6','Location','NorthEastOutside')
    legend boxoff
    ylabel('x')
    
    set(gcf,'Position',[100 300 1200 500])
end
end
