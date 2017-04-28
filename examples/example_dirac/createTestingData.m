function [ output_args ] = createTestingData( input_args )
    %CREATETESTINGDATA This function writes some data for continuous
    % integration tests
    
    !head -n -1 simulate_model_dirac.m > simulate_model_dirac_hdf.m
    !tail -n +2 writeSimulationData.template.m >> simulate_model_dirac_hdf.m
    
    global amiHDFfile
    global amiHDFprefix;
    amiHDFfile = fullfile(pwd, '../../tests/dirac/expectedResults.h5');
    
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

