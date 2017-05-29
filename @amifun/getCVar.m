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
        case 'JSparse'
            this.cvar =  'Jdata';
        case 'JB'
            this.cvar =  'JBdata';
        case 'JSparseB'
            this.cvar =  'JBdata';
        case 'w'
            this.cvar =  'udata->w';
        case 'dwdx'
            this.cvar =  'dwdx';
        case 'dwdp'
            this.cvar =  'udata->dwdp';
        case 'M'
            this.cvar =  'udata->M';
        case 'dfdx'
            this.cvar =  'udata->dfdx';
        case 'sz_tf'
            this.cvar =  'sz';
        otherwise
            this.cvar = this.funstr;
    end
end

