function [ this ] = getDeps(this, model)
    % getDeps populates the sensiflag for the requested function
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
            
        case 'dsigma_zdp'
            this.deps = {'sigma_z','p'};
            
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
            this.deps = {'deltax','deltaxdot','ddeltaxdx','ddeltaxdp','ddeltaxdt','dtaudp','xdot','sx','stau'};
        
        case 'deltaqB'
            this.deps = {'ddeltaxdp','xB'};
            
        case 'deltaxB'
            this.deps = {'deltax','dtaudp','xdot','xB'};
            
        case 'z'
            this.deps = {'x','k','p'};
            
        case 'dzdp'
            this.deps = {'z','p'};
            
        case 'dzdx'
            this.deps = {'z','x'};
            
        case 'dzdt'
            this.deps = {'z','x','xdot'};
            
        case 'sz'
            this.deps = {'dzdp','dzdx','dzdt','sx','dtaudp'};
            
        case 'sz_tf'
            this.deps = {'dzdp','dzdx','sx'};
            
        case 'dtaudp'
            this.deps = {'drootdp','drootdt'};
            
        case 'stau'
            this.deps = {'sroot','drootdt'};
            
        case 'sroot'
            this.deps = {'drootdp','drootdx','sx'};
 
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
        
        case 'sigma_z'
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
        
        case 'deltaxdot'
            this.deps = {};
            
        case 'Jy'
            this.deps = {'y','sigma_y','p'};
        case 'dJydx'
            this.deps = {'Jy','x'};
        case 'dJydsigma'
            this.deps = {'Jy','sigma_y'};
        case 'dJydp'
            this.deps = {'Jy','p','dJydsigma'};
        case 'sJy'
            this.deps = {'dJydp','dJydx','sy'};
            
        case 'Jz'
            this.deps = {'z','sigma_z','p'};
        case 'dJzdx'
            this.deps = {'Jz','x'};
        case 'dJzdsigma'
            this.deps = {'Jz','sigma_z'};
        case 'dJzdp'
            this.deps = {'Jz','p','dJzdsigma'};
        case 'sJz'
            this.deps = {'dJzdp','dJydx','sz'};


        otherwise
            error(['unknown function string: ' this.funstr ])
    end
end

