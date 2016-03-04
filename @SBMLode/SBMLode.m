classdef SBMLode < handle
    %SBMLMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        state@sym;
        observable@sym;
        parameter@sym;
        constant@sym;
        reaction@sym;
        compartment@sym;
        volume@sym;
        initState@sym;
        condition@sym;
        flux@sym;
        stochiometry@sym;
        xdot@sym;
        trigger@sym;
        bolus@sym;
        funmath@cell;
        funarg@cell;
        p@sym;
        k@sym;
    end
    
    methods
        function model = SBMLode(filename)
            model.importSBML(filename);
        end
        
        
        importSBML(this,filename)
        writeAMICI(this,modelname)
    end
    
end

