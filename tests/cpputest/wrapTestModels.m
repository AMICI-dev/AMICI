function [ output_args ] = wrapTestModels( input_args )
    %WRAPTESTMODELS Summary of this function goes here
    %   Detailed explanation goes here
    
    amiciPath = fileparts(mfilename('fullpath'));
    amiciPath = [amiciPath '/../..'];

    %% EXAMPLE STEADYSTATE

    cd([amiciPath '/examples/example_steadystate/']);

    [exdir,~,~]=fileparts(which('example_steadystate.m'));
    amiwrap('model_steadystate','model_steadystate_syms',exdir)
          
    %% EXAMPLE DIRAC
    cd([amiciPath '/examples/example_dirac/']);
    
    [exdir,~,~]=fileparts(which('example_dirac.m'));
    amiwrap('model_dirac','model_dirac_syms',exdir)
    
    cd(fileparts(mfilename('fullpath')));
    
end

