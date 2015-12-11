function [ this ] = getDeps(this, model)
    % getDeps writes the dependencies for the requested function
    %
    % Parameters:
    %  model: model definition object @type amimodel
    %
    % Return values:
    %  this: updated function definition object @type amifun
    
    switch(this.funstr)
        case 'xdot'
            if(strcmp(model.wtype,'iw'))
                this.deps = {'M','p','x','k','dx'};
            else
                this.deps = {'p','x','k'};
            end
            
        case 'dfdx'
            this.deps = {'xdot','x'};
            
        case 'J'
            if(strcmp(model.wtype,'iw'))
                this.deps = {'dfdx','M','x','xdot'};
            else
                this.deps = {'xdot','x'};
            end
            
        case 'JB'
            this.deps = {'J'};
            
        case 'dxdotdp'
            this.deps = {'xdot','p'};
            
        case 'sx0'
            this.deps = {'x0','p'};    
            
        case 'sdx0'
            this.deps = {'dx0','p'}; 
            
        case 'sxdot'
            if(strcmp(model.wtype,'iw'))
                this.deps = {'dfdx','M','dxdotdp','sdx','sx'};
            else
                this.deps = {'J','dxdotdp','sx'};
            end
            
        case 'dydx'
            this.deps = {'y','x'};
            
        case 'dydp'
            this.deps = {'y','p'};
          
        case 'sy'
            this.deps = {'dydp','dydx','sx'};
            
        case 'Jv'
            this.deps = {'J'};
            
        case 'JvB'
            this.deps = {'J'};
            
        case 'xBdot'
            if(strcmp(model.wtype,'iw'))
                this.deps = {'J','M','xB','dxB'};
            else
                this.deps = {'J','xB'};
            end
            
        case 'qBdot'
            this.deps = {'dxdotdp','xB'};
            
        case 'dsigma_ydp'
            this.deps = {'sigma_y','p'};
            
        case 'dsigma_tdp'
            this.deps = {'sigma_t','p'};
            
        case 'root'
            this.deps = {'x','k','p'};
            
        case 'drootdp'
            this.deps = {'root','p','drootdx','sx'};
            
        case 'drootdx'
            this.deps = {'root','x'};
            
        case 'drootdt'
            this.deps = {'root','drootdx','xdot'};
            
        case 'deltax'
            this.deps = {'x','k','p'};
            
        case 'ddeltaxdp'
            this.deps = {'deltax','p'};
            
        case 'ddeltaxdx'
            this.deps = {'deltax','x'};
            
        case 'ddeltaxdt'
            this.deps = {'deltax'};
            
        case 'deltasx'
            this.deps = {'deltax','dtaudp','xdot','sx','deltaxdot'};
        
        case 'deltaqB'
            this.deps = {'deltax','dtaudp','xdot','sx'};
            
        case 'deltaxB'
            this.deps = {'deltax','dtaudp','xdot','sx'};
            
        case 'z'
            this.deps = {'x','k','p'};
            
        case 'dzdp'
            this.deps = {'z','p'};
            
        case 'dzdx'
            this.deps = {'z','x'};
            
        case 'sz'
            this.deps = {'dzdp','dzdx','sx','dtaudp'};
            
        case 'x0'
            this.deps = {'p','k'};
            
        case 'JBand'
            this.deps = {'J'};
            
        case 'JBandB'
            this.deps = {'JB'};
            
        case 'JSparse'
            this.deps = {'J'};
            
        case 'JSparseB'
            this.deps = {'JB'};
            
        case 'y'
            this.deps = {'x','p','k'};
            
        case 'sigma_y'
            this.deps = {'p','k'};
        
        case 'sigma_t'
            this.deps = {'p','k'};
            
        case 'rhs'
            this.deps = {'xdot'};
            
        case 'dx0'
            this.deps = {'x','p','k'};
            
        case 'M'
            this.deps = {'x','p','k'};
            
        case 'x'
            this.deps = {};
            
        case 'dx'
            this.deps = {};
            
        case 'xB'
            this.deps = {};
            
        case 'dxB'
            this.deps = {};
            
        case 'k'
            this.deps = {};
            
        case 'p'
            this.deps = {};
            
        case 'sx'
            this.deps = {};
            
        case 'sdx'
            this.deps = {};
            
            
        otherwise
            error(['unknown function string: ' funstr ])
    end
end

