function [ output_args ] = createTestingData( input_args )
    %CREATETESTINGDATA This function writes some data for continuous
    % integration tests
    
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
    
end

