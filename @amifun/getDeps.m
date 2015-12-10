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
                deps = {'M','p','x','k','dx'};
            else
                deps = {'p','x','k'};
            end
            
        case 'dfdx'
            deps = {'xdot','x'};
            
        case 'J'
            if(strcmp(this.wtype,'iw'))
                deps = {'dfdx','M','x','xdot'};
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
                deps = {'dfdx','M','dxdotdp','sdx','sx'};
            else
                deps = {'J','dxdotdp','sx'};
            end
            
        case 'dydx'
            deps = {'y','x'};
            
        case 'dydp'
            deps = {'y','p'};
          
        case 'sy'
            deps = {'dydp','dydx','sx'};
            
        case 'drootdx'
            deps = {'rfun','x'};
            
        case 'drootdt'
            deps = {'rfun','drootdx','xdot'};
            
        case 'drootpdp'
            deps = {'rfun','p'};
            
        case 'drootdp'
            if(this.nxtrue > 0)
                deps = {'drootpdp','drootdx','dxdp','sx'};
            else
                deps = {'drootpdp','drootdx','sx'};
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
            deps = {'rfun','x'};
            
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
                deps = {'J','M','xB','dxB'};
            else
                deps = {'J','xB'};
            end
            
        case 'qBdot'
            deps = {'dxdotdp','xB'};
            
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
            deps = {'p','k'};
            
        case 'JBand'
            deps = {'J'};
            
        case 'JBandB'
            deps = {'JB'};
            
        case 'JSparse'
            deps = {'J'};
            
        case 'JSparseB'
            deps = {'JB'};
            
        case 'y'
            deps = {'x','p','k'};
            
        case 'rootval'
            deps = {'rfun'};
            
        case 'drvaldx'
            deps = {'drootdx'};
            
        case 'drvaldp'
            deps = {'drootdp'};
            
        case 'root'
            deps = {'rfun','rdisc'};
            
        case 'sigma_y'
            deps = {'p','k'};
        
        case 'sigma_t'
            deps = {'p','k'};
            
        case 'rhs'
            deps = {'xdot'};
            
        case 'dx0'
            deps = {'x','p','k'};
            
        case 'rfun'
            deps = {'p','k','x'};
            
        case 'M'
            deps = {'x','p','k'};
            
        case 'x'
            deps = {};
            
        case 'dx'
            deps = {};
            
        case 'xB'
            deps = {};
            
        case 'dxB'
            deps = {};
            
        case 'k'
            deps = {};
            
        case 'p'
            deps = {};
            
        case 'sx'
            deps = {};
            
        case 'sdx'
            deps = {};
            
            
        otherwise
            error(['unknown function string: ' funstr ])
    end
end

