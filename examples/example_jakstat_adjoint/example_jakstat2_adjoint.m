function example_jakstat2_adjoint()
    
    % compile the model
    [exdir,~,~]=fileparts(which('example_jakstat2_adjoint.m'));
    % amiwrap('model_jakstat2_adjoint','model_jakstat2_adjoint_syms',exdir,2)
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
    xi_rand = xi - 0.2;
    runs = 50;
    options.atol = 1e-13;
    options.rtol = 1e-13;
    
    % Get time for simulation
    tic;
    for i = 1 : runs    
        options.sensi = 0;
        sol0 = simulate_model_jakstat2_adjoint([],xi_rand,[],D,options);
    end
    t0 = toc;
    grad_sol0 = zeros(17,1);
    
    for j = 1:17
        epsvec = zeros(17,1);
        epsscal= 1e-5;
        epsvec(j) = epsscal;
        sol0_f = simulate_model_jakstat2_adjoint([],xi_rand + epsvec,[],D,options);
        sol0_b = simulate_model_jakstat2_adjoint([],xi_rand - epsvec,[],D,options);
        grad_sol0(j) = (sol0_f.llh - sol0_b.llh) / (2*epsscal);
    end
    
    % Get time for usual evaluation
    tic;
    for i = 1 : runs    
        options.sensi = 1;
        options.sensi_meth = 'adjoint';
        sol1 = simulate_model_jakstat2_adjoint([],xi_rand,[],D,options);
        grad_sol1 = sol1.sllh;
    end
    t1 = toc;
    
    % Get time for Finite Differences
    hvp = zeros(17,1);
    hvp_f = zeros(17,1);
    hvp_b = zeros(17,1);
    tic;
    for i = 1 : runs
        sol2 = simulate_model_jakstat2_adjoint([],xi_rand,[],D,options);
        v = sol2.sllh;
        delta = 1e-5;
        solp  = simulate_model_jakstat2_adjoint([],xi_rand + delta*v,[],D,options);
        solm  = simulate_model_jakstat2_adjoint([],xi_rand - delta*v,[],D,options);
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
        presol = simulate_model_jakstat2_adjoint([],xi_rand,[],D,options);
        v = presol.sllh;
        options.sensi = 2;
        sol  = simulate_model_jakstat2_adjoint([],xi_rand,[],D,options,v);
        hvpasa = hvpasa + sol.s2llh;
        grad_sol = sol.sllh;
    end
    t3 = toc;
    hvpasa = hvpasa / runs;
    hvptable = [hvp, hvp_f, hvp_b, hvpasa, abs(hvpasa - hvp) ./ hvp, abs(hvpasa - hvp) ./ hvpasa];
    gradtable = [grad_sol0, grad_sol1, grad_sol, abs(grad_sol0 - grad_sol1) ./ grad_sol0, ...
        abs(grad_sol0 - grad_sol1) ./ grad_sol1, abs(grad_sol - grad_sol1) ./ grad_sol1];  
    
    fprintf('Time elapsed for %i ODE solves (no sensis): %12.7f \n\n', runs, t0);
    fprintf('Time elapsed for %i gradient computations:  %12.7f \n\n', runs, t1);
    fprintf('Time elapsed for %i HVP computations (FD):  %12.7f \n\n', runs, t2);
    fprintf('Time elapsed for %i HVP computations (ASA): %12.7f \n\n', runs, t3);
    
    fprintf('| HVP from FD (c) | HVP from FD (f) | HVP from FD (b) |  HVP from ASA   | rel.Err.b.o.FD  | rel.Err.b.o.ASA |\n');
    fprintf('|=================|=================|=================|=================|=================|=================|\n');
    for i = 1 : 17
       fprintf('| %15.9f | %15.9f | %15.9f | %15.9f | %15.9f | %15.9f |\n', hvptable(i,:));
    end
    fprintf('|===========================================================================================================|\n');
    
    fprintf('|  Grad from FD   | Grad from ASAo1 | Grad from ASAo2 | rel.Err.b.o.FD  | rel.Err.b.o.ASA | r.E.ASAo2/ASAo1 |\n');
    fprintf('|=================|=================|=================|=================|=================|=================|\n');
    for i = 1 : 17
       fprintf('| %15.9f | %15.9f | %15.9f | %15.9f | %15.9f | %15.9f |\n', gradtable(i,:));
    end
    fprintf('|===========================================================================================================|\n');
    
    
end
