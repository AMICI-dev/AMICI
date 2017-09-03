function createTestingData()
% createTestingData runs simulation on all test models and writes results as hdf5. currently necessary for continuous integration
% to yield meaningful results
%
% Return values:
%  void

    
    oldwd = pwd;
    
    amiciPath = fileparts(mfilename('fullpath'));
    amiciPath = [amiciPath '/../..'];

    global amiHDFfile
    global amiHDFprefix;
    amiHDFfile = [amiciPath '/tests/cpputest/expectedResults.h5'];
    delete(amiHDFfile) 

    %% EXAMPLE STEADYSTATE

    cd([amiciPath '/examples/example_steadystate/']);

    [exdir,~,~]=fileparts(which('example_steadystate.m'));
    amiwrap('model_steadystate','model_steadystate_syms',exdir)
          
    !sed '$d' simulate_model_steadystate.m > simulate_model_steadystate_hdf.m
    !tail -n +2 ../../tests/cpputest/writeSimulationData.template.m >> simulate_model_steadystate_hdf.m
        
    % time vector
    t = linspace(0,100,50);
    p = [1;0.5;0.4;2;0.1];
    k = [0.1,0.4,0.7,1];
    
    options = amioption('sensi',0,...
        'maxsteps',1e4);
        
    amiHDFprefix = '/model_steadystate/nosensi/';
    sol = simulate_model_steadystate_hdf([t,inf],log10(p),k,[],options);

    amiHDFprefix = '/model_steadystate/sensiforward/';
    options.sensi = 1;
    sol = simulate_model_steadystate_hdf([t,inf],log10(p),k,[],options);
    
    amiHDFprefix = '/model_steadystate/sensiforwarddense/';
    options.sensi = 1;
    options.linsol = 1;
    sol = simulate_model_steadystate_hdf([t,inf],log10(p),k,[],options);
    
    amiHDFprefix = '/model_steadystate/nosensiSPBCG/';
    options.sensi = 0;
    options.linsol = 7;
    sol = simulate_model_steadystate_hdf([t,inf],log10(p),k,[],options);

    %% EXAMPLE DIRAC
    cd([amiciPath '/examples/example_dirac/']);
    
    [exdir,~,~]=fileparts(which('example_dirac.m'));
    amiwrap('model_dirac','model_dirac_syms',exdir)

    !sed '$d' simulate_model_dirac.m > simulate_model_dirac_hdf.m
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
    
    %% EXAMPLE JAKSTAT
    cd([amiciPath '/examples/example_jakstat_adjoint/']);
    
    [exdir,~,~]=fileparts(which('example_jakstat_adjoint.m'));
     amiwrap('model_jakstat_adjoint', 'model_jakstat_adjoint_syms', exdir, 1);

    !sed '$d' simulate_model_jakstat_adjoint.m > simulate_model_jakstat_adjoint_hdf.m
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

    amiHDFprefix = '/model_jakstat_adjoint/sensi2forward/';
    options.sensi = 2;
    options.sensi_meth = 'forward';
    simulate_model_jakstat_adjoint_hdf([],xi_rand,[],D,options);
    
    amiHDFprefix = '/model_jakstat_adjoint/sensi2adjoint/';
    options.sensi = 2;
    options.sensi_meth = 'adjoint';
    simulate_model_jakstat_adjoint_hdf([],xi_rand,[],D,options);
    
    amiHDFprefix = '/model_jakstat_adjoint/sensiforwardlogparam/';
    options.sensi = 1;
    options.pscale = 'log';
    options.sensi_meth = 'forward';
    simulate_model_jakstat_adjoint_hdf([],log(10.^xi_rand),[],D,options);
    
    amiHDFprefix = '/model_jakstat_adjoint/sensi2forwardlogparam/';
    options.sensi = 2;
    options.pscale = 'log';
    options.sensi_meth = 'forward';
    simulate_model_jakstat_adjoint_hdf([],xi_rand,[],D,options);
    
    %% EXAMPLE NEURON
    cd([amiciPath '/examples/example_neuron/']);
    
    [exdir,~,~]=fileparts(which('example_neuron.m'));
     amiwrap('model_neuron', 'model_neuron_syms', exdir, 1);

    !sed '$d' simulate_model_neuron.m > simulate_model_neuron_hdf.m
    !tail -n +2 ../../tests/cpputest/writeSimulationData.template.m >> simulate_model_neuron_hdf.m

    p = [0.02 0.3 65 0.9];

    k = [-60,10];
    
    t = linspace(0,100,100);

    options = amioption('sensi',0,...
        'maxsteps',1e5,...
        'nmaxevent',22);

    sol = simulate_model_neuron(t,log10(p),k,[],options);
    
    D = amidata(length(t),1,size(sol.z(sol.z<t(end)),2),size(sol.z(sol.z<t(end)),1),2);
    
    rng(0);
    D.Z = sol.z(sol.z<t(end));
    t = linspace(0,D.Z(end)-0.1,100);
    D.Z = D.Z + 0.5*randn(size(D.Z));
    D.Z(3) = NaN;
    D.Sigma_Z = 0.5*ones(size(D.Z));
    D.Z = D.Z + D.Sigma_Z.*randn(size(D.Z));
    
    D.t = t;
    D.condition = k;
    
    amiHDFprefix = '/model_neuron/nosensi/';
    options.sensi = 0;
    simulate_model_neuron_hdf([],log10(p),[],D,options);
    
    options.sensi = 1;
    amiHDFprefix = '/model_neuron/sensiforward/';
    options.sensi_meth = 'forward';
    simulate_model_neuron_hdf([],log10(p),[],D,options);

    amiHDFprefix = '/model_neuron/sensi2forward/';
    options.sensi = 2;
    options.sensi_meth = 'forward';
    simulate_model_neuron_hdf([],log10(p),[],D,options);
    
    
    %% EXAMPLE EVENTS
    cd([amiciPath '/examples/example_events/']);
    
    [exdir,~,~]=fileparts(which('example_events.m'));
     amiwrap('model_events', 'model_events_syms', exdir);

    !sed '$d' simulate_model_events.m > simulate_model_events_hdf.m
    !tail -n +2 ../../tests/cpputest/writeSimulationData.template.m >> simulate_model_events_hdf.m

    p = [0.5;2;0.5;0.5];

    k = [4,8,10,4];
    
    t = linspace(0,10,20);

    options = amioption('sensi',0,...
        'maxsteps',1e4);
    D = amidata(length(t),1,2,2,4);
    
    D.t = t;
    D.condition = k;
    
    amiHDFprefix = '/model_events/nosensi/';
    options.sensi = 0;
    simulate_model_events_hdf([],log10(p),k,D,options);
    
    amiHDFprefix = '/model_events/sensiforward/';
    options.sensi = 1;
    options.sensi_meth = 'forward';
    simulate_model_events_hdf([],log10(p),[],D,options);
    
    %% EXAMPLE NESTED_EVENTS
    cd([amiciPath '/examples/example_nested_events/']);
    
    [exdir,~,~]=fileparts(which('example_nested_events.m'));
    amiwrap('model_nested_events', 'model_nested_events_syms', exdir);
    
    !sed '$d' simulate_model_nested_events.m > simulate_model_nested_events_hdf.m
    !tail -n +2 ../../tests/cpputest/writeSimulationData.template.m >> simulate_model_nested_events_hdf.m
    
    t = linspace(0,20,100);
    p = [ 0.1
        1000
        2
        8e-1
        1.6   ];
    k = [];
    
    options = amioption('sensi',0,...
        'maxsteps',1e4,...
        'rtol',1e-14,...
        'ato',1e-12,...
        'nmaxevent', 2);
    D = amidata(length(t),1,0,2,0);
    
    amiHDFprefix = '/model_nested_events/nosensi/';
    options.sensi = 0;
    simulate_model_nested_events_hdf(t,log10(p),k,D,options);
    
    amiHDFprefix = '/model_nested_events/sensiforward/';
    options.sensi = 1;
    options.sensi_meth = 'forward';
    simulate_model_nested_events_hdf(t,log10(p),k,D,options);

    %%
    chdir(oldwd)
end

