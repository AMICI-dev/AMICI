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
    runs = 1;
    options.atol = 1e-12;
    options.rtol = 1e-12;
    
    % Get time for simulation
    tic;
    for i = 1 : runs    
        options.sensi = 0;
        sol0 = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options);
    end
    t0 = toc;
    
    % Get time for usual evaluation
    tic;
    for i = 1 : runs    
        options.sensi = 1;
        options.sensi_meth = 'adjoint';
        sol1 = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options);
    end
    t1 = toc;
    
    % Get time for Finite Differences
    hvp = zeros(17,1);
    hvp_f = zeros(17,1);
    hvp_b = zeros(17,1);
    tic;
    for i = 1 : runs
        sol2 = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options);
        v = sol2.sllh;
        delta = 1e-5;
        solp  = simulate_model_jakstat_adjoint_hvp([],xi_rand + delta*v,[],D,options);
        solm  = simulate_model_jakstat_adjoint_hvp([],xi_rand - delta*v,[],D,options);
        hvp = hvp + (solp.sllh - solm.sllh) / (2*delta);
        hvp_f = hvp_f + (solp.sllh - sol2.sllh) / (delta);
        hvp_b = hvp_b + (sol2.sllh - solm.sllh) / (delta);
    end
    t2 = toc;
    hvp = hvp / runs;
    hvp_f = hvp_f / runs;
    hvp_b = hvp_b / runs;
    
    % Get time for Second order adjoints
    hvpasa = zeros(17,1);
    tic;
    for i = 1 : runs
        options.sensi = 1;
        options.sensi_meth = 'adjoint';
        presol = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options);
        v = presol.sllh;
        options.sensi = 2;
        sol  = simulate_model_jakstat_adjoint_hvp([],xi_rand,[],D,options,v);
        hvpasa = hvpasa + sol.s2llh;
    end
    t3 = toc;
    
%     % Printed output
%     hvpasa = hvpasa / runs;
%     hvptable = [hvp, hvp_f, hvp_b, hvpasa, abs(hvpasa - hvp) ./ hvp, abs(hvpasa - hvp) ./ hvpasa];
%     fprintf('Time elapsed for %i ODE solves (no sensis): %12.7f \n\n', runs, t0);
%     fprintf('Time elapsed for %i gradient computations:  %12.7f \n\n', runs, t1);
%     fprintf('Time elapsed for %i HVP computations (FD):  %12.7f \n\n', runs, t2);
%     fprintf('Time elapsed for %i HVP computations (ASA): %12.7f \n\n', runs, t3);
%     
%     fprintf('| HVP from FD (c) | HVP from FD (f) | HVP from FD (b) |  HVP from ASA   | rel.Err.b.o.FD  | rel.Err.b.o.ASA |\n');
%     fprintf('|=================|=================|=================|=================|=================|=================|\n');
%     for i = 1 : 17
%        fprintf('| %15.9f | %15.9f | %15.9f | %15.9f | %15.9f | %15.9f |\n', hvptable(i,:));
%     end
%     fprintf('|===========================================================================================================|\n');
    
    figure();
    
    subplot(1,2,1);
    hold on;
    min_x = min(hvpasa);
    max_x = max(hvpasa);
    plot([min_x, max_x], [min_x, max_x],'k:');
    plot(hvpasa(1:14), hvp_f(1:14), 'ro');
    plot(hvpasa(1:14), hvp_b(1:14), 'go');
    plot(hvpasa(1:14), hvp(1:14), 'bo');
    title('Hessian vector product with gradient');
    xlabel('Second order adjoints sensitivities');
    ylabel('Finite differences');
    legend('', 'Forward FD', 'Backward FD', 'Centered FD');
    box on;
    hold off;
    
    subplot(1,2,2);
    hold on;
    plot(1, t0, 'mo');
    plot(2, t1, 'ro');
    plot(3, t2, 'bo');
    plot(4, t3, 'ko');
    title(['Runtime for ' num2str(runs) 'evaluations (in sec)']);
    legend('ODE Integration', 'Gradient computation (ASA)', 'HVP from FD via 1st order ASA', 'HVP via 2nd order ASA');
    box on;
    hold off;
end
