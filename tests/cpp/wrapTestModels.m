function wrapTestModels()
% wrapTestModels calls amiwrap on all test models. currently necessary for continuous integrations
% to yield meaningful results
%
% Return values:
%  void
    
    amiciPath = fileparts(mfilename('fullpath'));
    amiciPath = [amiciPath '/../../matlab'];
    
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
    
    %% EXAMPLE EVENTS
    cd([amiciPath '/examples/example_events/']);
    
    try
        [exdir,~,~]=fileparts(which('example_events.m'));
        amiwrap('model_events', 'model_events_syms', exdir);
    catch err
        disp(err.message)
        cd(fileparts(mfilename('fullpath')));
    end
    
    %% EXAMPLE NESTED EVENTS
    cd([amiciPath '/examples/example_nested_events/']);
    
    try
        [exdir,~,~]=fileparts(which('example_nested_events.m'));
        amiwrap('model_nested_events', 'model_nested_events_syms', exdir);
    catch err
        disp(err.message)
        cd(fileparts(mfilename('fullpath')));
    end

    cd(fileparts(mfilename('fullpath')));
    
    %% EXAMPLE ROBERTSON
    cd([amiciPath '/examples/example_robertson/']);
    
    try
        [exdir,~,~]=fileparts(which('example_robertson.m'));
        amiwrap('model_robertson', 'model_robertson_syms', exdir);
    catch err
        disp(err.message)
        cd(fileparts(mfilename('fullpath')));
    end

    cd(fileparts(mfilename('fullpath')));
    
    %% EXAMPLE CALVETTI
    cd([amiciPath '/examples/example_calvetti/']);
    
    try
        [exdir,~,~]=fileparts(which('example_calvetti.m'));
        amiwrap('model_calvetti', 'model_calvetti_syms', exdir);
    catch err
        disp(err.message)
        cd(fileparts(mfilename('fullpath')));
    end

    cd(fileparts(mfilename('fullpath')));
    
end

