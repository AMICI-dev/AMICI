%
% @file amidata
% @brief definition of optsym class
%
classdef optsym<sym
    %OPTSYM is an auxiliary class to gain access to the private symbolic
    %property 's' which is necessary to be able to call symobj::optimize 
    %on it

    properties
    end
    
    methods
        function obj=optsym(symbol)
            %optsym converts the symbolic object into a optsym object
            obj=obj@sym(symbol);
            %
            % Parameters:
            %  symbol: symbolic object @type sym
            %
            % Return values:
        end
        function out=getoptimized(obj)    
            %optsym calls symobj::optimize on the optsym object
            %
            % Parameters:
            %
            % Return values:
            %  optimized symbolic object @type sym
            out = mupadmex('symobj::optimize',obj.s);
        end
    end
    
end
