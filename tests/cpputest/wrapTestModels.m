function wrapTestModels()
% wrapTestModels calls amiwrap on all test models. currently necessary for continuous integrations
% to yield meaningful results
    
    amiciPath = fileparts(mfilename('fullpath'));
    amiciPath = [amiciPath '/../..'];
    
    %% EXAMPLE STEADYSTATE
    
    cd([amiciPath '/examples/example_steadystate/']);
    
    try
        [exdir,~,~]=fileparts(which('example_steadystate.m'));
        amiwrap('model_steadystate','model_steadystate_syms',exdir);
    catch err
        disp(err.message)
        cd(fileparts(mfilename('fullpath')));
    end
    
    %% EXAMPLE DIRAC
    cd([amiciPath '/examples/example_dirac/']);
    
    try
        [exdir,~,~]=fileparts(which('example_dirac.m'));
        amiwrap('model_dirac','model_dirac_syms',exdir);
    catch err
        disp(err.message)
        cd(fileparts(mfilename('fullpath')));
    end
    
    %% EXAMPLE JAKSTAT
    cd([amiciPath '/examples/example_jakstat_adjoint/']);
    
    try
        [exdir,~,~]=fileparts(which('example_jakstat_adjoint.m'));
        amiwrap('model_jakstat_adjoint', 'model_jakstat_adjoint_syms', exdir, 1);
    catch err
        disp(err.message)
        cd(fileparts(mfilename('fullpath')));
    end

    %% EXAMPLE NEURON
    cd([amiciPath '/examples/example_neuron/']);
    
    try
        [exdir,~,~]=fileparts(which('example_neuron.m'));
        amiwrap('model_neuron', 'model_neuron_syms', exdir, 1);
    catch err
        disp(err.message)
        cd(fileparts(mfilename('fullpath')));
    end

    cd(fileparts(mfilename('fullpath')));
    
end

