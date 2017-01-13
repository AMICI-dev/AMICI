function example_jakstat_adjoint()
    
    % compile the model
    [exdir,~,~]=fileparts(which('example_jakstat_adjoint.m'));
    amiwrap('model_jakstat','model_jakstat_adjoint_syms',exdir)
    
    num = xlsread(fullfile(exdir,'pnas_data_original.xls'));
    
    D.t = num(:,1);
    D.condition= [1.4,0.45];
    D.Y = num(:,[2,4,6]);
    D.Sigma_Y = NaN(size(D.Y));
    D = amidata(D);
    
    xi =  [0.60
        3
        -0.95
        -0.0075
        0
        -2.8
        -0.26
        -0.075
        -0.41
        -5
        -0.74
        -0.64
        -0.11
        0.027
        -0.5
        0
        -0.5];
    
    options.sensi = 0;
    sol = simulate_model_jakstat([],xi,[],D,options);
    
    if(usejava('jvm'))
    figure
    for iy = 1:3
        subplot(2,2,iy)
        plot(D.t,D.Y(:,iy),'rx')
        hold on
        plot(sol.t,sol.y(:,iy),'.-')
        xlim([0,60])
        xlabel('t')
        switch(iy)
            case 1
                ylabel('pStat')
            case 2
                ylabel('tStat')
            case 3
                ylabel('pEpoR')
        end
        ylim([0,1.2])
    end
    set(gcf,'Position',[100 300 1200 500])
    end
    
    % generate new
    xi_rand = xi + 0.1;
    options.sensi = 1;
    options.sensi_meth = 'adjoint';
    sol = simulate_model_jakstat([],xi_rand,[],D,options);
    options.sensi_meth = 'forward';
    solf = simulate_model_jakstat([],xi_rand,[],D,options);
    
    options.sensi = 0;
    eps = 1e-4;
    fd_grad = NaN(length(xi),1);
    for ip = 1:length(xi)
        xip = xi_rand;
        xip(ip) = xip(ip) + eps;
        psol = simulate_model_jakstat([],xip,[],D,options);
        fd_grad(ip) = (psol.llh-sol.llh)/eps;
    end
    
    if(usejava('jvm'))
    figure
    bar([abs((sol.sllh-fd_grad)./sol.sllh),abs((sol.sllh-solf.sllh)./sol.sllh)])
    hold on
    set(gca,'YScale','log')
    ylim([1e-12,1e0])
    box on
    hold on
%     plot([1e-2,1e2],[1e-2,1e2],'k:')
    xlabel('parameter index')
    ylabel('relative difference to adjoint sensitivities')
    legend('finite differences','forward sensitivities')
    set(gcf,'Position',[100 300 1200 500])
    end
    
    
    drawnow
    
end
