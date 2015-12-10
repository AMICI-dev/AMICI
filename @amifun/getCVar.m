function [ cvar ] = getFunCVar(this, funstr )
    % getDeps returns the c variable for the requested function.
    %
    % Parameters:
    %  funstr: function or expression for which the dependencies should be returned  @type string
    %
    % Return values:
    %  cvar: cell array of dependencies @type string
    
    switch(funstr)
        case 'J'
            cvar =  'Jdata';
        case 'rootval'
            cvar =  'rootval';
        case 'JSparse'
            cvar =  'Jdata';
        case 'JB'
            cvar =  'JBdata';
        case 'JSparseB'
            cvar =  'JBdata';
        case 'drvaldp'
            cvar =  'drvaldp';
        case 'drvaldx'
            cvar =  'drvaldx';
        otherwise
            cvar = funstr;
    end
end

