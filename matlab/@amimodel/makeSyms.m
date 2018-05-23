function makeSyms( this )
% makeSyms extracts symbolic definition from the user provided model
% and checks them for consistency
%
% Parameters:
%
% Return values:
%  void



% check whether sym is properly defined
if(~isfield(this.sym,'xdot'))
    if(isfield(this.sym,'f'))
        try
            this.sym.xdot = betterSym(this.sym.f(:));
            this.sym = rmfield(this.sym,'f');
        catch
            error('Could not transform model.sym.f into a symbolic variable, please check the definition!')
        end
    else
        error('field xdot/f is missing in model definition')
    end
else
    if(isfield(this.sym,'f'))
        if(~isequaln(this.sym.f,this.sym.xdot))
            error('Model this contains conflicting definitions sym.f and sym.xdot of DE right hand side');
        end
    else
        try
            if(~isa(this.sym.xdot(:),'sym'))
                this.sym.xdot = betterSym(this.sym.xdot(:));
            else
                this.sym.xdot = this.sym.xdot(:);
            end
        catch
            error('Could not transform model.sym.xdot into a symbolic variable, please check the definition!')
        end
    end
end

if(~isfield(this.sym,'x'))
    error('Model this is missing the definition of the state vector x (.sym.x)!')
else
    try
        this.sym.x = sym(this.sym.x(:));
    catch
        error('Could not transform model.sym.x into a symbolic variable, please check the definition!')
    end
end
if(numel(this.sym.x)~=numel(this.sym.xdot))
    error('Size of model.sym.x and model.sym.xdot does not agree.')
end

if(~isfield(this.sym,'p'))
    error('Model this is missing the definition of the parameter vector p (.sym.p)!')
else
    try
        this.sym.p = sym(this.sym.p(:));
    catch
        error('Could not transform model.sym.y into a symbolic variable, please check the definition!')
    end
end
if(~isfield(this.sym,'x0'))
    error('Model this is missing the definition of the vector of initial conditions x0 (.sym.x0)!')
else
    try
        if(~isa(this.sym.x0(:),'sym'))
            this.sym.x0 = betterSym(this.sym.x0(:));
        else
            this.sym.x0 = this.sym.x0(:);
        end
    catch
        error('Could not transform model.sym.x0 into a symbolic variable, please check the definition!')
    end
end
if(numel(this.sym.x)~=numel(this.sym.x0))
    error('Size of model.sym.x and model.sym.x0 does not agree.')
end
if(any(ismember(symvar(this.sym.x0),this.sym.x)))
    error('initial states x0 must not contain state variables x');
end

if(~isfield(this.sym,'y'))
    error('Model this is missing the definition of the vector of observables y (.sym.y)!')
else
    try
        this.sym.y = sym(this.sym.y(:));
    catch
        error('Could not transform model.sym.y into a symbolic variable, please check the definition!')
    end
end

% complete optional fields
if(~isfield(this.sym,'u'))
    this.sym.u = sym(zeros(0,0));
end
if(~isfield(this.sym,'k'))
    this.sym.k = sym(zeros(0,0));
else
    try
        this.sym.k = sym(this.sym.k(:));
    catch
        error('Could not transform model.sym.k into a symbolic variable, please check the definition!')
    end
    
end

if(isfield(this.sym,'root'))
    error('The definition of events via a root function is deprecated and no longer supported. Please update the model definition syntax!')
end
if(~isfield(this.sym,'sigma_y'))
    this.sym.sigma_y = sym(ones(length(this.sym.y),1));
end
if(numel(this.sym.sigma_y) == 1)
    this.sym.sigma_y = this.sym.sigma_y*sym(ones(length(this.sym.y),1));
end
if(numel(this.sym.sigma_y)~=length(this.sym.y))
    error('Size of model.sym.y and model.sym.sigma_y does not agree.')
end

if(any(ismember(this.sym.k,this.sym.p)))
    error(['Invalid Model: ' char(this.sym.k(find(ismember(this.sym.k,this.sym.p),1))) ' is contained in both p and k!'])
end

if(~isfield(this.sym,'Jy'))
    this.sym.Jy = sym(zeros(numel(this.sym.y),1));
    for iy = 1:length(this.sym.y)
        this.sym.Jy(iy) = betterSym(['0.5*log(2*pi*sigma_y_' num2str(iy-1) '^2) + 0.5*((y_' num2str(iy-1) '-my_' num2str(iy-1) ')/sigma_y_' num2str(iy-1) ')^2']);
    end
end

symvars = symvar(this.sym.xdot);
for ivar = 1:length(symvars)
    if(isequaln(symvars(ivar),sym('E')))
        error('The symbolic entities named ''E'' are currently not supported due to restrictions of the symbolic math toolbox!');
    end
end
%     svaridx = not(ismember(symvars,[this.sym.p;this.sym.k;this.sym.x;sym('t')]));
%     if(any(svaridx))
%         error(['The symbolic variable ' char(symvars(find(svaridx,1))) ' is used in the differential equation right hand side but was not specified as parameter/state/constant!']);
%     end
end

