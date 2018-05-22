function example_robertson()
%%
% COMPILATION

[exdir,~,~]=fileparts(which('example_robertson.m'));
% compile the model
amiwrap('model_robertson','model_robertson_syms',exdir)

%%
% SIMULATION

% time vector
t = [0 4*logspace(-6,6)];
p = [0.04;1e4;3e7];
k = [0.9];

options = amioption('sensi',0,...
    'maxsteps',1e4,...
    'atol',1e-6,...
    'rtol',1e-4);
% load mex into memory
[~] = which('simulate_model_robertson'); % fix for inaccessability problems
sol = simulate_model_robertson(t,log10(p),k,[],options);

tic
sol = simulate_model_robertson(t,log10(p),k,[],options);
disp(['Time elapsed with cvodes: ' num2str(toc) ])

%%
% ODE15S

dae_system = @(t,y,p,k) [-p(1)*y(1) + p(2)*y(2).*y(3)
   p(1)*y(1) - p(2)*y(2).*y(3) - p(3)*y(2).^2
   y(1) + y(2) + y(3) - 1 ];
y0 = [k(1); 0; 0];
M = [1 0 0
   0 1 0
   0 0 0];
options_ode15s = odeset('Mass',M,'RelTol',options.rtol,'AbsTol',options.atol,'MaxStep',options.maxsteps);

tic
[~, X_ode15s] = ode15s(@(t,x) dae_system(t,x,p,k),t,y0,options_ode15s);
disp(['Time elapsed with ode15s: ' num2str(toc) ])

%%
% PLOTTING
if(usejava('jvm'))
    figure
    c_x = get(gca,'ColorOrder');
    subplot(2,2,1)
    for ix = 1:size(sol.x,2)
        plot(t,sol.y(:,ix),'.-','Color',c_x(ix,:))
        hold on
        if(ix==2)
            scale = 1e4;
        else
            scale = 1;
        end
        plot(t,scale*X_ode15s(:,ix),'d','Color',c_x(ix,:))
    end
    legend('x1','x1_{ode15s}','x2','x2_{ode15s}','x3','x3_{ode15s}','Location','NorthEastOutside')
    legend boxoff
    xlabel('time t')
    set(gca,'XScale','log')
    ylabel('x')
    box on
    subplot(2,2,2)
    plot(t,abs(sol.x-X_ode15s),'--')
    set(gca,'YScale','log')
    legend('error x1','error x2','error x3','Location','NorthEastOutside')
    legend boxoff
    ylabel('x')
    
    set(gcf,'Position',[100 300 1200 500])
end

%%
% FORWARD SENSITIVITY ANALYSIS

options.sensi = 1;

sol = simulate_model_robertson(t,log10(p),k,[],options);

%%
% FINITE DIFFERENCES

eps = 1e-2;
xi = log10(p);
for ip = 1:length(p);
    xip = xi;
    xip(ip) = xip(ip) + eps;
    solp = simulate_model_robertson(t,xip,k,[],options);
    sy_fd(:,:,ip) = (solp.y - sol.y)/eps;
end

%%
% PLOTTING
if(usejava('jvm'))
    figure
    for ip = 1:length(p)
        subplot(length(p),2,ip*2-1)
        hold on
        for ix = 1:size(sol.x,2)
            plot(t,sol.sy(:,ix,ip),'.-','Color',c_x(ix,:))
            plot(t,sy_fd(:,ix,ip),'d','Color',c_x(ix,:))
        end
        legend('sy1','sy1_{fd}','sy2','sy2_{fd}','sy3','sy3_{fd}','Location','NorthEastOutside')
        legend boxoff
        title(['state sensitivity for p' num2str(ip)])
        xlabel('time t')
        set(gca,'XScale','log')
        ylabel('sx')
        box on
        
        subplot(length(p),2,ip*2)
        plot(t,abs(sol.sy(:,:,ip)-sy_fd(:,:,ip)),'--')
        legend('error sy1','error sy2','error sy3','Location','NorthEastOutside')
        legend boxoff
        title(['state sensitivity for p' num2str(ip)])
        xlabel('time t')
        set(gca,'XScale','log')
        ylabel('error')
        set(gca,'YScale','log')
        box on
    end
    set(gcf,'Position',[100 300 1200 500])
    
    
    drawnow
end
end