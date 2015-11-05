function [ deps ] = getFunDeps(this, funstr )
    % getDeps returns the dependencies for the requested function/expression.
    %
    % Parameters:
    %  funstr: function or expression for which the dependencies should be returned  @type string
    %
    % Return values:
    %  deps: cell array of dependencies @type cell
    
    switch(funstr)
        case 'xdot'
            if(strcmp(this.wtype,'iw'))
                deps = {'M','rhs'};
            else
                deps = {'xdot'};
            end
            
        case 'dfdx'
            deps = {'xdot','x'};
            
        case 'J'
            if(strcmp(this.wtype,'iw'))
                deps = {'dfdx','M'};
            else
                deps = {'xdot','x'};
            end
            
        case 'JB'
            deps = {'J'};
            
        case 'dxdotdp'
            deps = {'xdot','p'};
            
        case 'sx0'
            deps = {'x0','p'};    
            
        case 'sdx0'
            deps = {'dx0','p'}; 
            
        case 'sxdot'
            if(strcmp(this.wtype,'iw'))
                deps = {'dfdx','M','dxdotdp'};
            else
                deps = {'J','dxdotdp'};
            end
            
        case 'dydx'
            deps = {'y','x'};
            
        case 'dydp'
            deps = {'y','p'};
          
        case 'sy'
            deps = {'dydp','dydx'};
            
        case 'drootdx'
            deps = {'root','x'};
            
        case 'drootdt'
            deps = {'root','drootdx','xdot'};
            
        case 'drootpdp'
            deps = {'root','p'};
            
        case 'drootdp'
            if(this.nxtrue > 0)
                deps = {'drootpdp','drootdx','dxdp'};
            else
                deps = {'drootpdp','drootdx'};
            end
            
        case 'sroot'
            deps = {'drootdt','drootdp'};
            
        case 'srootval'
            deps = {'drootdp'};
            
        case 's2root'
            deps = {'drootdt','drootdp','sroot','ddrootdtdt','ddrootdtdp','s2rootval'};
            
        case 's2rootval'
            deps = {'ddrootdpdp'};
            
        case 'dxdp'
            deps = {'x'};
            
        case 'ddrootdxdx'
            deps = {'root','x'};
            
        case 'ddrootdxdt'
            deps = {'drootdx','xdot','x','ddrootdxdx'};
           
        case 'ddrootdtdt'
            deps = {'drootdt','ddrootdxdt','xdot'};
        
        case 'ddxdtdp'
            deps = {'dxdotdp','dxdotdx','dxdp'};
            
        case 'dxdotdx'
            deps = {'xdot','x'};
            
        case 'ddrootdpdx'
            deps = {'drootdp','x'};
        
        case 'ddrootdpdt'
            deps = {'drootdp','ddrootdpdx','xdot',};
            
        case 'ddrootdtdp'
            deps = {'drootdt','p','ddrootdxdt','dxdp'};
            
        case 'drootdx_ddxdpdp'
            deps = {'drootdx','p'};
            
        case 'ddrootdpdp'
            deps = {'drootdp','p','ddrootdpdx','dxdp','drootdx_ddxdpdp'};
            
        case 'dtdp'
            deps = {'drootdt','drootpdp'};
            
        case 'dtdx'
            deps = {'drootdt','drootdx'};
            
        case 'Jv'
            deps = {'J'};
            
        case 'JvB'
            deps = {'J'};
            
        case 'xBdot'
            if(strcmp(this.wtype,'iw'))
                deps = {'J','M'}; %% TBD
            else
                deps = {'J'};
            end
            
        case 'qBdot'
            deps = {'dxdotdp'};
            
        case 'dsigma_ydp'
            deps = {'sigma_y','p'};
            
        case 'dsigma_tdp'
            deps = {'sigma_t','p'};
            
        case 'rdisc'
            deps = {'xdot','x'};
            
        case 'deltadisc'
            deps = {'rdisc','xdot','x'};
            
        case 'ideltadisc'
            deps = {'deltadisc','rdisc','xdot','qBdot','x'};
            
        case 'bdeltadisc'
            deps = {'deltadisc','rdisc','xdot','xBdot','x'};
            
        case 'sdeltadisc'
            deps = {'deltadisc','rdisc','xdot','sxdot','x'};
            
        case 'x0'
            deps = {'x0'};
            
        case 'JBand'
            deps = {'J'};
            
        case 'JBandB'
            deps = {'JB'};
            
        case 'JSparse'
            deps = {'J'};
            
        case 'JSparseB'
            deps = {'JB'};
            
        case 'y'
            deps = {'y'};
            
        case 'rootval'
            deps = {'rfun'};
            
        case 'drvaldx'
            deps = {'drootdx'};
            
        case 'drvaldp'
            deps = {'drootdp'};
            
        case 'root'
            deps = {'rfun','rdisc'};
            
        case 'sigma_y'
            deps = {'sigma_y'};
        
        case 'sigma_t'
            deps = {'sigma_t'};
            
        case 'rhs'
            deps = {'xdot'};
            
        case 'dx0'
            deps = {'dx0'};
            
        case 'rfun'
            deps = {'rfun'};
            
        otherwise
            error(['unknown function string: ' funstr ])
    end
end

