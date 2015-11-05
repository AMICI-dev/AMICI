function [ svar ] = getFunSVar(this, funstr )
    % getDeps returns the symbolic variable for the requested function.
    %
    % Parameters:
    %  funstr: function or expression for which the dependencies should be returned  @type string
    %
    % Return values:
    %  cvar: cell array of dependencies @type string
   
    switch(funstr)
        case 'rootval'
            svar =  'root';
        case 'drvaldp'
            svar =  'drootpdp';
        case 'drvaldx'
            svar =  'drootdx';
        otherwise
            svar = funstr;
    end
end

