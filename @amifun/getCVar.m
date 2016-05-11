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
            this.cvar =  'w_tmp';
        case 'dwdx'
            this.cvar =  'dwdx_tmp';
        case 'dwdp'
            this.cvar =  'dwdp_tmp';
        case 'M'
            this.cvar =  'M_tmp';
        case 'dfdx'
            this.cvar =  'dfdx_tmp';
        otherwise
            this.cvar = this.funstr;
    end
end

