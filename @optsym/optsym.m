classdef optsym<sym

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
