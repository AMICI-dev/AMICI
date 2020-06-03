function example_steadystate
    %%
    % COMPILATION
    
    [exdir,~,~]=fileparts(which('example_steadystate.m'));
    % compile the model
    amiwrap('model_steadystate','model_steadystate_syms',exdir)
    
    %%
    % SIMULATION
    
    % time vector
    t = linspace(0,100,50);
    p = [1;0.5;0.4;2;0.1];
    k = [0.1,0.4,0.7,1];
    
    options = amioption(...
        'sensi', 0, ...
        'maxsteps', 1e4 ...
        );
    
    % load mex into memory
    simulate_model_steadystate(t,log10(p),k,[],options);
    
    tic;
    sol = simulate_model_steadystate([t, inf],log10(p),k,[],options);
    display(['Time elapsed with cvodes: ' num2str(toc) ' seconds']);
    
    %%
    % ODE15S
    
    ode_system = @(t,x,p,k) [-2*p(1)*x(1)^2 - p(2)*x(1)*x(2) + 2*p(3)*x(2) + p(4)*x(3) + p(5);
        + p(1)*x(1)^2 - p(2)*x(1)*x(2) - p(3)*x(2) + p(4)*x(3);
        + p(2)*x(1)*x(2) - p(4)*x(3) - k(4)*x(3)];
    options_ode15s = odeset('RelTol',options.rtol,'AbsTol',options.atol,'MaxStep',options.maxsteps);
    
    tic;
    [~, X_ode15s] = ode15s(@(t,x) ode_system(t,x,p,k),t,k(1:3),options_ode15s);
    disp(['Time elapsed with ode15s: ' num2str(toc) ' seconds'])
    
    %%
    % PLOTTING
    
    if(usejava('jvm'))
        figure('Name', 'Example SteadyState');
        c_x = get(gca,'ColorOrder');
        subplot(2,1,1);
        for ix = 1:size(sol.x,2)
            plot(t,sol.x(1:end-1,ix),'.-','Color',c_x(ix,:));
            hold on;
            plot(t,X_ode15s(:,ix),'d','Color',c_x(ix,:));
            plot([t(1), t(end)],sol.x(end,ix)*[1, 1],'--','Color',c_x(ix,:));
        end
        legend('x1','x1_{ode15s}','x1_{ss, Newton}','x2','x2_{ode15s}','x2_{ss, Newton}','x3','x3_{ode15s}','x3_{ss, Newton}','Location','NorthEastOutside');
        legend boxoff;
        xlabel('time t');
        ylabel('x');
        box on;
        subplot(2,1,2);
        plot(t,abs(sol.x(1:end-1,:)-X_ode15s),'--');
        set(gca,'YScale','log');
        legend('error x1','error x2','error x3','Location','NorthEastOutside');
        legend boxoff;
        set(gcf,'Position',[100 300 1200 500]);
    end
    
    %%
    % FORWARD SENSITIVITY ANALYSIS
    
    options.sensi = 1;
    options.sens_ind = [3,1,2,4];
    sol = simulate_model_steadystate([t, inf],log10(p),k,[],options);
    
    %%
    % FINITE DIFFERENCES
    
    eps = 1e-3;
    
    xi = log10(p);
    sx_ffd = zeros(length(t)+1, 3, length(p));
    sx_bfd = zeros(length(t)+1, 3, length(p));
    sx_cfd = zeros(length(t)+1, 3, length(p));
    for ip = 1:4;
        xip = xi;
        xim = xi;
        xip(ip) = xip(ip) + eps;
        xim(ip) = xim(ip) - eps;
        solp = simulate_model_steadystate([t, inf],xip,k,[],options);
        solm = simulate_model_steadystate([t, inf],xim,k,[],options);
        sx_ffd(:,:,ip) = (solp.x - sol.x) / eps;
        sx_bfd(:,:,ip) = (sol.x - solm.x) / eps;
        sx_cfd(:,:,ip) = (solp.x - solm.x) / (2*eps);
    end
    
    %%
    % PLOTTING
    if(usejava('jvm'))
        % Sensitivities for time series
        figure('Name', 'Example SteadyState');
        for ip = 1:4
            subplot(4,2,ip*2-1);
            hold on;
            for ix = 1:size(sol.x,2)
                plot(t,sol.sx(1:end-1,ix,ip),'.-','Color',c_x(ix,:));
                plot(t,sx_cfd(1:end-1,ix,options.sens_ind(ip)),'d','Color',c_x(ix,:));
            end
            legend('x1','x1_{fd}','x2','x2_{fd}','x3','x3_{fd}','Location','NorthEastOutside');
            legend boxoff;
            title(['state sensitivity for p' num2str(options.sens_ind(ip))]);
            xlabel('time t');
            ylabel('x');
            box on;

            subplot(4,2,ip*2);
            plot(t,abs(sol.sx(1:end-1,:,ip)-sx_cfd(1:end-1,:,options.sens_ind(ip))),'--');
            legend('error x1','error x2','error x3','Location','NorthEastOutside');
            legend boxoff;
            title(['error of state sensitivity for p' num2str(options.sens_ind(ip))]);
            xlabel('time t');
            ylabel('error');
            set(gca,'YScale','log');
            box on;
        end
        set(gcf,'Position',[100 300 1200 500]);
        
        sxss = squeeze(sol.sx(length(t),:,:));
        sxss_fd = squeeze(sx_cfd(length(t),:,options.sens_ind));
        
        % Sensitivities for steady state
        figure('Name', 'Example SteadyState');
        subplot(1,2,1);
        hold on;
        plot([-1,1], [-1,1], 'k:');
        type = {'o','x','d','*'};
        for ip = 1:4
            for ix = 1:size(sol.x,2)
                plot(sxss(ix,ip), sxss_fd(ix,ip), type{ip}, 'Color', c_x(ix,:));
            end
        end
        legend('SA = FD',...
            's^{x1}_{p3}','s^{x2}_{p3}','s^{x3}_{p3}',...
            's^{x1}_{p1}','s^{x2}_{p1}','s^{x3}_{p1}',...
            's^{x1}_{p2}','s^{x2}_{p2}','s^{x3}_{p2}',...
            's^{x1}_{p4}','s^{x2}_{p3}','s^{x3}_{p4}',...
            'Location','NorthEastOutside');
        legend boxoff;
        title('Steady state sensitivies');
        xlabel('Steady state sensitivities');
        ylabel('finite differences');
        box on;
        
        
        subplot(1,2,2);
        hold on;
        for ip = 1:4
            for ix = 1:size(sol.x,2)
                ind = 3*(ip-1)+ix;
                plot(ind, abs(sxss(ind)-sxss_fd(ind)), type{ip}, 'Color', c_x(ix,:));
            end
        end
        title('Error of steady state sensitivies');
        set(gca,'XTick',1:12);
        set(gca,'XTickLabel',...
            {'s^{x1}_{p3}','s^{x2}_{p3}','s^{x3}_{p3}',...
            's^{x1}_{p1}','s^{x2}_{p1}','s^{x3}_{p1}',...
            's^{x1}_{p2}','s^{x2}_{p2}','s^{x3}_{p2}',...
            's^{x1}_{p4}','s^{x2}_{p3}','s^{x3}_{p4}'});
        xlabel('Sensitivity index');
        xlim([0.5,12.5])
        ylabel('Error');
        box on;
        set(gca,'YScale','log');
        set(gcf,'Position',[100 300 1200 500]);
    end
    
    %%
    % XDOT FOR DIFFERENT TIME POINTS
    
    t = [10,25,100,250,1000];
    options.sensi = 0;
    ssxdot = NaN(length(t), size(sol.x, 2));
    for it = 1:length(t)
        tt = [0,t(it)];
        solss = simulate_model_steadystate(tt,log10(p),k,[],options);
        ssxdot(it,:) = solss.diagnosis.xdot;
    end

    % Compute steady state wihtout integration before
    sol = simulate_model_steadystate(inf,log10(p),k,[],options);
    

    % Test recapturing in the case of Newton solver failing
    options.newton_maxsteps = 4;
    options.maxsteps = 300;
    sol_newton_fail = simulate_model_steadystate(inf,log10(p),k,[],options);
    
    %%
    % PLOTTING
    if(usejava('jvm'))
        figure('Name', 'Example SteadyState');
        subplot(1,3,[1 2]);
        hold on;
        for ix = 1:size(sol.x,2)
            plot(t,abs(ssxdot(:,ix)),'o-','Color',c_x(ix,:));
            plot([10 1000], abs(sol_newton_fail.diagnosis.xdot(ix))*[1, 1],'--','Color',c_x(ix,:));
        end
        plot(t,sqrt(sum(ssxdot.^2,2)), 'ko-');
        plot([10 1000],sqrt(sum(sol_newton_fail.diagnosis.xdot.^2,2))*[1, 1], 'k--');
        legend('dx_1/dt','dx_1/dt, ss','dx_2/dt','dx_2/dt, ss',...
            'dx_3/dt','dx_3/dt, ss','||dx/dt||','||dx/dt||, ss',...
            'Location','NorthEastOutside');
        legend boxoff;
        title('ODE right hand side over time');
        xlabel('time t');
        ylabel('xdot');
        box on;
        set(gca,'YScale','log');
        set(gca,'XScale','log');
        
        subplot(1,3,3);
        hold on;
        bar(sol_newton_fail.diagnosis.posteq_numsteps([1, 3]));
        legend boxoff;
        title('Number of Newton steps');
        xlabel('Solver run');
        ylabel('Newton steps');
        xlim([0.5,2.5]);
        ylim([-1 5]);
        a = gca();
        a.Children.BarWidth = 0.6;
        box on;
        
        set(gcf,'Position',[100 300 1200 500]);
    end
    
end


