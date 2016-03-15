%
% @file amioption
% @brief definition of amioption class
%
classdef SBMLode < handle
    %SBMLODE carries all information about the differential equation
    %defined by a SBML model definition file. This class acts as an
    %interface between SBML files and \ref amimodel
    
    properties
        % vector of non-constant and non-boundary states @type symbolic
        state = sym([]);
        % vector of guessed observables @type symbolic
        observable = sym([]);
        % vector of guessed observable names @type symbolic
        observable_name = sym([]);
        % vector of SBML parameters @type symbolic
        param = sym([]);
        % vector of amimodel parameters @type symbolic
        parameter = sym([]);
        % vector of constant states @type symbolic
        constant = sym([]);
        % vector of compartments @type symbolic
        compartment = sym([]);
        % vector of compartment volumes @type symbolic
        volume = sym([]);
        % vector of initial values for non-constant and non-boundary states @type symbolic
        initState = sym([]);
        % vector of boundary condition states which are not constant @type symbolic
        condition = sym([]);
        % vector of reaction fluxes @type symbolic
        flux = sym([]);
        % matrix of reaction stochiometries @type symbolic
        stochiometry = sym([]);
        % right hand side of differential equation for states @type symbolic
        xdot = sym([]);
        % vector of trigger functions for events @type symbolic
        trigger = sym([]);
        % matrix of event bolus @type symbolic
        bolus = sym([]);
        % cell array containing the function definition @type cell
        funmath = {};
        % cell array containing the function arguments @type cell
        funarg = {};
        % expression for time in the model @type char
        time_symbol = char([]);
        % nominal parameter value @type sym
        pnom = double([]);
        % nominal condition value @type sym
        knom = double([]);
    end
    
    methods
        function model = SBMLode(filename)
            model.importSBML(filename);
        end
        
        
        importSBML(this,filename)
        writeAMICI(this,modelname)
    end
    
end

