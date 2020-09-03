function [ this ] = getCVar(this)
    % getCVar populates the cvar property
    %
    % Parameters:
    %
    % Return values:
    %  this: updated function definition object @type amifun
    
    switch(this.funstr)
        case 'JSparse'
            this.cvar =  'var_JSparse';
        case 'M'
            this.cvar =  'M';
        case 'dfdx'
            this.cvar =  'dfdx';
        otherwise
            this.cvar = ['var_' this.funstr];
    end
end

