function [ output_args ] = wrapTestModels( input_args )
    %WRAPTESTMODELS Summary of this function goes here
    %   Detailed explanation goes here
    
    amiciPath = fileparts(mfilename('fullpath'));
    amiciPath = [amiciPath '/../..'];
    
    %% EXAMPLE STEADYSTATE
    
    cd([amiciPath '/examples/example_steadystate/']);
    
    try
        [exdir,~,~]=fileparts(which('example_steadystate.m'));
        amiwrap('model_steadystate','model_steadystate_syms',exdir);
    catch
        cd(fileparts(mfilename('fullpath')));
    end
    
    %% EXAMPLE DIRAC
    cd([amiciPath '/examples/example_dirac/']);
    
    try
        [exdir,~,~]=fileparts(which('example_dirac.m'));
        amiwrap('model_dirac','model_dirac_syms',exdir);
    catch
        cd(fileparts(mfilename('fullpath')));
    end
    
    cd(fileparts(mfilename('fullpath')));
    
end

