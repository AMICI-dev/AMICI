function [ output_args ] = createTestingData( input_args )
    %CREATETESTINGDATA This function writes some data for continuous
    % integration tests
    
    oldwd = pwd;
    
    amiciPath = fileparts(mfilename('fullpath'));
    amiciPath = [amiciPath '/../..'];

    global amiHDFfile
    global amiHDFprefix;
    amiHDFfile = [amiciPath '/tests/cpputest/expectedResults.h5'];

    %% EXAMPLE STEADYSTATE

    cd([amiciPath '/examples/example_steadystate/']);

    [exdir,~,~]=fileparts(which('example_steadystate.m'));
    amiwrap('model_steadystate','model_steadystate_syms',exdir)
          
    !head -n -1 simulate_model_steadystate.m > simulate_model_steadystate_hdf.m
    !tail -n +2 ../../tests/cpputest/writeSimulationData.template.m >> simulate_model_steadystate_hdf.m
        
    % time vector
    t = linspace(0,300,20);
    p = [1;0.5;0.4;2;0.1];
    k = [0.1,0.4,0.7,1];
    
    options = amioption('sensi',0,...
        'maxsteps',1e4);
        
    amiHDFprefix = '/model_steadystate/nosensi/';
    sol = simulate_model_steadystate_hdf(t,log10(p),k,[],options);

    amiHDFprefix = '/model_steadystate/sensiforward/';
    options.sensi = 1;
    sol = simulate_model_steadystate_hdf(t,log10(p),k,[],options);

    %% EXAMPLE DIRAC
    cd([amiciPath '/examples/example_dirac/']);
    
    [exdir,~,~]=fileparts(which('example_dirac.m'));
    amiwrap('model_dirac','model_dirac_syms',exdir)

    !head -n -1 simulate_model_dirac.m > simulate_model_dirac_hdf.m
    !tail -n +2 ../../tests/cpputest/writeSimulationData.template.m >> simulate_model_dirac_hdf.m
       
    % time vector
    t = linspace(0,3,1001);
    p = [1;0.5;2;3];
    k = [];    
    options = amioption('sensi',0,...
        'maxsteps',1e4);
        
    amiHDFprefix = '/model_dirac/nosensi/';
    sol = simulate_model_dirac_hdf(t,log10(p),k,[],options);
    
    amiHDFprefix = '/model_dirac/sensiforward/';
    options.sensi = 1;
    sol = simulate_model_dirac_hdf(t,log10(p),k,[],options);
    
    %% EXAMPLE JAKSTAT ADJOINT
    cd([amiciPath '/examples/example_jakstat_adjoint/']);
    
    [exdir,~,~]=fileparts(which('example_jakstat_adjoint.m'));
     amiwrap('model_jakstat_adjoint', 'model_jakstat_adjoint_syms', exdir, 1);

    !head -n -1 simulate_model_jakstat_adjoint.m > simulate_model_jakstat_adjoint_hdf.m
    !tail -n +2 ../../tests/cpputest/writeSimulationData.template.m >> simulate_model_jakstat_adjoint_hdf.m

    num = xlsread(fullfile(exdir,'pnas_data_original.xls'));
    
    D.t = num(:,1);
    D.condition = [1.4,0.45];
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
    
    amiHDFprefix = '/model_jakstat_adjoint/nosensi/';
    options.sensi = 0;
    simulate_model_jakstat_adjoint_hdf([],xi,[],D,options);
    
    xi_rand = xi + 0.1;
    options.sensi = 1;
    amiHDFprefix = '/model_jakstat_adjoint/sensiadjoint/';
    options.sensi_meth = 'adjoint';
    simulate_model_jakstat_adjoint_hdf([],xi_rand,[],D,options);
    
    amiHDFprefix = '/model_jakstat_adjoint/sensiforward/';
    options.sensi_meth = 'forward';
    simulate_model_jakstat_adjoint_hdf([],xi_rand,[],D,options);

    amiHDFprefix = '/model_jakstat_adjoint/sensi2adjoint/';
    amiwrap('model_jakstat_adjoint', 'model_jakstat_adjoint_syms', exdir, 1);
    options.sensi = 2;
    options.sensi_meth = 'forward';
    simulate_model_jakstat_adjoint_hdf([],xi_rand,[],D,options);

    %%
    chdir(oldwd)
end

