%
% @file optsym
% @brief this is a placeholder class to get access to the protected sym.s
% property which is necessary to apply the mupad command symobj::optimize
%
classdef optsym<sym
    %OPTSYM is a placeholder class to get access to the protected sym.s 

    properties
    end
    
    methods
        function obj=optsym(symbol)
            %optsym creates an optsym object from the symbolic expression
            %symbol
            %
            % Parameters:
            %  symbol: symbolic expression @type sym
            % 
            % Return values:
            %  obj: transformed symbolic expression
           obj=obj@sym(symbol); 
        end
        function out=getoptimized(obj)    
            % getoptimized generates an optimized formulation of the
            % expression of the symbolic expression obj. this optimized
            % formulation identifies frequently employed expressions and
            % subsitutes them by placeholder variables to avoid repeated
            % evaluation
            % 
            % Return values:
            %  out: optimized expression
            out = mupadmex('symobj::optimize',obj.s);
        end
    end
    
end
