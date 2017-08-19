function example_jakstat_adjoint()
    
    % compile the model
    [exdir,~,~]=fileparts(which('example_jakstat_adjoint.m'));
    amiwrap('model_jakstat_adjoint','model_jakstat_adjoint_syms',exdir,1)
    
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
    sol = simulate_model_jakstat_adjoint([],xi,[],D,options);
    
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
    options.sensi = 2;
    options.sensi_meth = 'adjoint';
    sol = simulate_model_jakstat_adjoint([],xi_rand,[],D,options);
    options.sensi_meth = 'forward';
    solf = simulate_model_jakstat_adjoint([],xi_rand,[],D,options);
    
    options.sensi = 1;
    eps = 1e-4;
    fd_grad = NaN(length(xi),1);
    for ip = 1:length(xi)
        xip = xi_rand;
        xip(ip) = xip(ip) + eps;
        psol = simulate_model_jakstat_adjoint([],xip,[],D,options);
        fd_grad(ip) = (psol.llh-sol.llh)/eps;
        fd_hess(:,ip) = (psol.sllh-sol.sllh)/eps;
    end
    
    if(usejava('jvm'))
    figure
    subplot(1,2,1)
    plot(abs(solf.sllh),abs(fd_grad),'rx')
    hold on
    plot(abs(solf.sllh),abs(sol.sllh),'bo')
    hold on
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    ylim([1e-2,1e2])
    xlim([1e-2,1e2])
    plot([1e-2,1e2],[1e-2,1e2],'k:')
    legend('finite differences','adjoint sensitivities','Location','best')
    box on
    axis square
    xlabel('absolute value forward sensitivity gradient entries')
    ylabel('absolute value gradient entries')
    
    subplot(1,2,2)
    plot(abs(solf.s2llh(:)),abs(fd_hess(:)),'rx')
    hold on
    plot(abs(solf.s2llh(:)),abs(sol.s2llh(:)),'bo')
    hold on
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    ylim([1e-5,1e3])
    xlim([1e-5,1e3])
    plot([1e-5,1e3],[1e-5,1e3],'k:')
    legend('finite differences','adjoint sensitivities','Location','best')
    box on
    axis square
    xlabel('absolute value forward sensitivity hessian entries')
    ylabel('absolute value hessian entries')
    
    set(gcf,'Position',[100 300 1200 500])
    end
    
    
    drawnow
    
end
