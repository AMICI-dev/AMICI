function [ this ] = getSVar(this)
    % getDeps returns the symbolic variable for the requested function.
    %
    % Parameters:
    %  funstr: function or expression for which the dependencies should be returned  @type string
    %
    % Return values:
    %  cvar: cell array of dependencies @type string
   
    switch(this.funstr)
        case 'rootval'
            this.svar =  'root';
        case 'drvaldp'
            this.svar =  'drootpdp';
        case 'drvaldx'
            this.svar =  'drootdx';
        otherwise
            this.svar = this.funstr;
    end
end

