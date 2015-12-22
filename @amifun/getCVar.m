function [ this ] = getCVar(this)
    % getCVar populates the cvar property
    %
    % Parameters:
    %
    % Return values:
    %  this: updated function definition object @type amifun
    
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
        otherwise
            this.cvar = this.funstr;
    end
end

