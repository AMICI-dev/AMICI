classdef SBMLode < handle
    %SBMLMODEL provides an intermediate container between the SBML
    %definition and an amimodel object
    
    properties
        % states
        state = sym.empty();
        % observables
        observable = sym.empty();
        % names of observables
        observable_name = sym.empty();
        % parameter names
        param = sym.empty();
        % parameter expressions 
        parameter = sym.empty();
        % constants
        constant = sym.empty();
        % reactions
        reaction = sym.empty();
        % compartments
        compartment = sym.empty();
        % compartment volumes
        volume = sym.empty();
        % condition volumes
        kvolume = sym.empty();
        % initial condition of states
        initState = sym.empty();
        % condition
        condition = sym.empty();
        % reaction fluxes
        flux = sym.empty();
        % reaction stochiometry
        stochiometry = sym.empty();
        % right hand side of reconstructed differential equation
        xdot = sym.empty();
        % event triggers
        trigger = sym.empty();
        % event boli
        bolus = sym.empty();
        % mathematical experessions for function
        funmath = cell.empty();
        % string function signature
        funarg = cell.empty();
        % symbol of time
        time_symbol = char.empty();
        % nominal parameters
        pnom = double.empty();
        % nominal conditions
        knom = double.empty();
    end
    
    methods
        function model = SBMLode(filename)
            % SBMLode extracts information from an SBML definition and
            % stores it in a symbolic format
            %
            % Parameters:
            %  filename: target name of the model (excluding the suffix
            %  .xml/.sbml)
            %
            % Return values:
            %
            model.importSBML(filename);
            model.checkODE();
        end
        
        importSBML(this,filename)

        checkODE(this)
        
        writeAMICI(this,modelname)
        
    end
    
end

