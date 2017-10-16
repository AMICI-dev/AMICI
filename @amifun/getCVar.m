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
        case 'Jv'
            this.cvar =  'Jv_tmp';
        case 'JvB'
            this.cvar =  'JvB_tmp';
        case 'JSparse'
            this.cvar =  'Jdata';
        case 'JB'
            this.cvar =  'JBdata';
        case 'JSparseB'
            this.cvar =  'JBdata';
        case 'M'
            this.cvar =  'tdata->M';
        case 'dfdx'
            this.cvar =  'tdata->dfdx';
        otherwise
            this.cvar = ['var_' this.funstr];
    end
end

