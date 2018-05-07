%
% @file optsym.m
% @brief wrapper class for sym to make symob::optimize accessible
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
            %
            % Parameters:
            %  symbol: symbolic object @type sym
            %
            % Return values:
            %  obj: cast symbolic object @type optsym
            obj=obj@sym(symbol);
        end
        function out=getoptimized(obj)    
            %getoptimized calls symobj::optimize on the optsym object
            %
            % Parameters:
            %
            % Return values:
            %  out: optimized symbolic object @type sym
            out = mupadmex('symobj::optimize',obj.s);
        end
    end
    
end
