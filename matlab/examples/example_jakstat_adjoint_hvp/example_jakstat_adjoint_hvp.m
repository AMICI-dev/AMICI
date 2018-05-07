function example_jakstat_adjoint_hvp()
    
    % compile the model
    [exdir,~,~]=fileparts(which('example_jakstat_adjoint_hvp.m'));
    amiwrap('model_jakstat_adjoint_hvp','model_jakstat_adjoint_hvp_syms',exdir,2)
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
    
    
    % generate new
    xi_rand = xi - 0.1;
    options.atol = 1e-12;
    options.rtol = 1e-12;
    
    % Get time for simulation
    tic;
    options.sensi = 0;
    sol0 = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options);
    t0 = toc;
    
    % Get time for usual evaluation
    tic;
    options.sensi = 1;
    options.sensi_meth = 'adjoint';
    sol1 = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options);
    t1 = toc;
    
    % Get time for Finite Differences
    hvp = zeros(17,1);
    hvp_f = zeros(17,1);
    hvp_b = zeros(17,1);
    tic;

    sol2 = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options);
    v = sol2.sllh;
    delta = 1e-5;
    solp  = simulate_model_jakstat_adjoint_hvp([],xi_rand + delta*v,[],D,options);
    solm  = simulate_model_jakstat_adjoint_hvp([],xi_rand - delta*v,[],D,options);
    hvp = hvp + (solp.sllh - solm.sllh) / (2*delta);
    hvp_f = hvp_f + (solp.sllh - sol2.sllh) / (delta);
    hvp_b = hvp_b + (sol2.sllh - solm.sllh) / (delta);
    t2 = toc;
    
    % Get time for Second order adjoints
    hvpasa = zeros(17,1);
    tic;
    options.sensi = 1;
    options.sensi_meth = 'adjoint';
    presol = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options);
    v = presol.sllh;
    options.sensi = 2;
    sol  = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options,v);
    hvpasa = hvpasa + sol.s2llh';
    t3 = toc;
    options.sensi_meth = 'forward';
    solf  = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options,v);
    options.sensi_meth = 'adjoint';

if(usejava('jvm'))
    figure();
    
    subplot(1,2,1);
    bar([abs((sol.s2llh-hvp)./sol.s2llh),abs((sol.s2llh-hvp_f)./sol.s2llh),abs((sol.s2llh-hvp_b)./sol.s2llh),abs((sol.s2llh-solf.s2llh)./sol.s2llh)])
    hold on
    set(gca,'YScale','log')
    ylim([1e-16,1e0])
    box on
    hold on
    %     plot([1e-2,1e2],[1e-2,1e2],'k:')
    xlabel('parameter index')
    ylabel('relative difference to adjoint sensitivities')
    legend('FD_{central}','FD_{forward}','FD_{backward}','forward sensitivities','Orientation','horizontal')
    legend boxoff
    set(gcf,'Position',[100 300 1200 500])
    
    subplot(1,2,2);
    hold on;
    bar([t0,t1,t2,t3]);
    xlabel('runtime [s]')
    set(gca,'XTick',1:4,'XTickLabel',{'ODE Integration', 'Gradient computation (ASA)', 'HVP from FD via 1st order ASA', 'HVP via 2nd order ASA'},'XTickLabelRotation',20);
    
    box on;
    hold off;
end
end
