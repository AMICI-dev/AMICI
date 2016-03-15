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
           obj=obj@sym(symbol); 
        end
        function out=getoptimized(obj)          
            out = mupadmex('symobj::optimize',obj.s);
        end
    end
    
end
