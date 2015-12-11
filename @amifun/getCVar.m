function [ this ] = getCVar(this)
    % getDeps returns the c variable for the requested function.
    %
    % Parameters:
    %  funstr: function or expression for which the dependencies should be returned  @type string
    %
    % Return values:
    %  cvar: cell array of dependencies @type string
    
    switch(this.funstr)
        case 'J'
            this.cvar =  'Jdata';
        case 'rootval'
            this.cvar =  'rootval';
        case 'JSparse'
            this.cvar =  'Jdata';
        case 'JB'
            this.cvar =  'JBdata';
        case 'JSparseB'
            this.cvar =  'JBdata';
        case 'drvaldp'
            this.cvar =  'drvaldp';
        case 'drvaldx'
            this.cvar =  'drvaldx';
        otherwise
            this.cvar = this.funstr;
    end
end

