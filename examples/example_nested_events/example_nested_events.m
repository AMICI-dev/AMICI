function example_events()
%%
% COMPILATION

[exdir,~,~]=fileparts(which('example_nested_events.m'));
% compile the model
amiwrap('model_nested_events','model_nested_events_syms',exdir)

%%
% SIMULATION

% time vector
t = linspace(0,20,100);
p = [ 0.1
    1000
    2
    8e-1
    1.6   ];
k = [];

options = amioption('sensi',0,...
    'maxsteps',1e4,...
    'rtol',1e-12,...
    'ato',1e-10,...
    'nmaxevent', 2);
D = amidata(length(t),1,0,2,0);
% load mex into memory
[~] = which('simulate_model_events'); % fix for inaccessability problems
sol = simulate_model_nested_events(t,log10(p),k,D,options);

tic
sol = simulate_model_nested_events(t,log10(p),k,D,options);
disp(['Time elapsed with cvodes: ' num2str(toc) ])

%%
% ODE15S

sig = 1e-2;
delta_num = @(tau) exp(-1/2*(tau/sig).^2)/(sqrt(2*pi)*sig);

ode_system = @(t,x,p,k) +p(2)*delta_num(t-p(3)) + x*p(4)*heaviside(x-1) - x*p(5);

options_ode15s = odeset('RelTol',options.rtol,'AbsTol',options.atol,'MaxStep',options.maxsteps);

tic
[~, X_ode15s] = ode15s(@(t,x) ode_system(t,x,p,k),t,p(1),options_ode15s);
disp(['Time elapsed with ode15s: ' num2str(toc) ])

%%
% PLOTTING

if(usejava('jvm'))
    figure
    c_x = get(gca,'ColorOrder');
    subplot(1,2,1)
    for ix = 1:size(sol.x,2)
        plot(t,sol.x(:,ix),'.-','Color',c_x(ix,:))
        hold on
        plot(t,X_ode15s(:,ix),'d','Color',c_x(ix,:))
    end
    legend('x1','x1_{ode15s}','Location','NorthEastOutside')
    legend boxoff
    xlabel('time t')
    ylabel('x')
    set(gca,'YScale','log')
    box on
    subplot(1,2,2)
    plot(t,abs(sol.x-X_ode15s),'--')
    set(gca,'YScale','log')
    legend('error x1','Location','NorthEastOutside')
    legend boxoff
    ylabel('x')
    
    set(gcf,'Position',[100 300 1200 300])
end


%%
% FORWARD SENSITIVITY ANALYSIS

options.sensi = 1;

sol = simulate_model_nested_events(t,log10(p),k,D,options);

%%
% FINITE DIFFERENCES

eps = 1e-8;
xi = log10(p);
for ip = 1:5;
    xip = xi;
    xip(ip) = xip(ip) + eps;
    solp = simulate_model_nested_events(t,xip,k,D,options);
    sx_fd(:,:,ip) = (solp.x - sol.x)/eps;
end

%%
% PLOTTING
if(usejava('jvm'))
    figure
    for ip = 1:5
        subplot(5,2,ip*2-1)
        hold on
        for ix = 1:size(sol.x,2)
            plot(t,sol.sx(:,ix,ip),'.-','Color',c_x(ix,:))
            plot(t,sx_fd(:,ix,ip),'d','Color',c_x(ix,:))
        end
        legend('sx1','sx1_{fd}','Location','NorthEastOutside')
        legend boxoff
        title(['state sensitivity for p' num2str(ip)])
        xlabel('time t')
        ylabel('sx')
        box on
        
        subplot(5,2,ip*2)
        plot(t,abs(sol.sx(:,:,ip)-sx_fd(:,:,ip)),'--')
        legend('error sx1','Location','NorthEastOutside')
        legend boxoff
        title(['state sensitivity for p' num2str(ip)])
        xlabel('time t')
        ylabel('error')
        set(gca,'YScale','log')
        box on
    end
    set(gcf,'Position',[100 300 1200 500])
    
    drawnow
end
end