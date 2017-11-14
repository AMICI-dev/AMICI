function [ this ] = getCVar(this)
    % getCVar populates the cvar property
    %
    % Parameters:
    %
    % Return values:
    %  this: updated function definition object @type amifun
    
    switch(this.funstr)
        case 'J'
            this.cvar =  'J';
        case 'JSparse'
            this.cvar =  'var_J';
        case 'Jv'
            this.cvar =  'var_Jv';
        case 'JvB'
            this.cvar =  'var_JvB';
        case 'JB'
            this.cvar =  'JB';
        case 'JSparseB'
            this.cvar =  'var_JB';
        case 'M'
            this.cvar =  'M';
        case 'dfdx'
            this.cvar =  'dfdx';
        otherwise
            this.cvar = ['var_' this.funstr];
    end
end

