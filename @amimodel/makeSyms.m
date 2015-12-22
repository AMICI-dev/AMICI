function [ this ] = makeSyms( this )
    % makeSyms extracts symbolic definition from the user provided model
    % and checks them for consistency
    %
    % Return values:
    %  this: updated model definition object @type amimodel
    
    % check whether sym is properly defined
    if(~isfield(this.sym,'xdot'))
        if(isfield(this.sym,'f'))
            this.sym.xdot = this.sym.f(:);
            this.sym = rmfield(this.sym,'f');
        else
            error('field xdot/f is missing in model definition')
        end
    else
        if(isfield(this.sym,'f'))
            if(~isequaln(this.sym.f,this.sym.xdot))
                error('Model this contains conflicting definitions sym.f and sym.xdot of DE right hand side');
            end
        else
            this.sym.xdot = this.sym.xdot(:);
        end
    end
    
    if(~isfield(this.sym,'x'))
        error('Model this is missing the definition of the state vector x (.sym.x)!')
    else
        this.sym.x = this.sym.x(:);
    end
    
    if(~isfield(this.sym,'p'))
        error('Model this is missing the definition of the parameter vector p (.sym.p)!')
    else
        this.sym.p = this.sym.p(:);
    end
    if(~isfield(this.sym,'x0'))
        error('Model this is missing the definition of the vector of initial conditions x0 (.sym.x0)!')
    else
        this.sym.x0 = this.sym.x0(:);
    end
    if(~isfield(this.sym,'y'))
        error('Model this is missing the definition of the vector of observables y (.sym.y)!')
    else
        this.sym.y = this.sym.y(:);
    end
    if(~all([size(this.sym.x,2)==size(this.sym.xdot,2),size(this.sym.xdot,2)==size(this.sym.x0,2)]))
        error('Sizes of x0, xdot and x do not agree!')
    end
    
    % complete optional fields
    if(~isfield(this.sym,'u'))
        this.sym.u = sym.empty(0,0);
    end
    if(~isfield(this.sym,'k'))
        this.sym.k = sym.empty(0,0);
    end
        
    if(isfield(this.sym,'root'))
        error('The definition of events via a root function is deprecated and no longer supported. Please update the model definition syntax!')
    end
    if(~isfield(this.sym,'sigma_y'))
        this.sym.sigma_y = sym.ones(size(this.sym.y));
    end
    if(numel(this.sym.sigma_y) == 1)
        this.sym.sigma_y = this.sym.sigma_y*sym.ones(size(this.sym.y));
    end
    
    if(any(ismember(this.sym.k,this.sym.p)))
        error(['Invalid Model: ' char(this.sym.k(find(ismember(this.sym.k,this.sym.p),1))) ' is contained in both p and k!'])
    end
    
    if(~isfield(this.sym,'Jy'))
        for iy = 1:length(this.sym.y)
            this.sym.Jy(iy) = sym(['0.5*log(2*pi*sdy_' num2str(iy-1) '^2) + 0.5*((y_' num2str(iy-1) '-my_' num2str(iy-1) ')/sdy_' num2str(iy-1) ')^2']);
        end
    end
    
end

