function this = parseModel(this)
    % parseModel parses the this and computes all necessary symbolic expressions.
    %
    % Return values:
    %  this: updated model definition object @type amimodel
    
    % check whether sym is properly defined
    if(~isfield(this.sym,'x'))
        error('Model this is missing the definition of the state vector x (.sym.x)!')
    end
    if(~isfield(this.sym,'xdot') && ~isfield(this.sym,'f'))
        error('Model this is missing the definition of the right hand side (.sym.xdot) or (.sym.f)!')
    end
    if(isfield(this.sym,'f'))
        if(isfield(this.sym,'xdot'))
            if(~isequaln(this.sym.f,this.sym.xdot))
                error('Model this contains conflicting definitions sym.f and sym.xdot of DE right hand side');
            end
        else
            this.sym.xdot = this.sym.f;
        end
    end
    
    if(~isfield(this.sym,'p'))
        error('Model this is missing the definition of the parameter vector p (.sym.p)!')
    end
    if(~isfield(this.sym,'x0'))
        error('Model this is missing the definition of the vector of initial conditions x0 (.sym.x0)!')
    end
    if(~isfield(this.sym,'y'))
        error('Model this is missing the definition of the vector of observables y (.sym.y)!')
    end
    if(size(this.sym.x,1)<size(this.sym.x,2))
        this.sym.x = transpose(this.sym.x);
    end
    if(size(this.sym.xdot,1)<size(this.sym.xdot,2))
        this.sym.xdot = transpose(this.sym.xdot);
    end
    
    if(size(this.sym.x0,1)<size(this.sym.x0,2))
        this.sym.x0 = transpose(this.sym.x0);
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
    if(~isfield(this.sym,'root'))
        this.sym.root = sym.empty(0,1);
    end
    if(~isfield(this.sym,'sigma_y'))
        this.sym.sigma_y = sym.ones(size(this.sym.y));
    end
    if(numel(this.sym.sigma_y) == 1)
        this.sym.sigma_y = this.sym.sigma_y*sym.ones(size(this.sym.y));
    end
    if(~isfield(this.sym,'sigma_t'))
        this.sym.sigma_t = sym.ones(size(this.sym.y));
    end
    if(numel(this.sym.sigma_t) == 1)
        this.sym.sigma_t = this.sym.sigma_t*sym.ones(size(this.sym.y));
    end
    
    if(any(ismember(this.sym.k,this.sym.p)))
        error(['Invalid Model: ' char(this.sym.k(find(ismember(this.sym.k,this.sym.p),1))) ' is contained in both p and k!'])
    end
    
    % check whether we have a DAE or ODE
    if(isfield(this.sym,'M'))
        this.wtype = 'iw'; % DAE
    else
        this.wtype = 'cw'; % ODE
    end

    % load old hashes
    [this,HTable] = this.loadOldHashes();
    
    nx = length(this.sym.x);
    np = length(this.sym.p);
    nk = length(this.sym.k);
    ny = length(this.sym.y);
    %remove zero-roots
    ir = 1;
    while ir <= length(this.sym.root)
        if(isequaln(this.sym.root(ir),0))
            this.sym.root(ir) = [];
        else
            ir = ir + 1;
        end
    end
    nr = length(this.sym.root);
    
    this.nx = nx;
    this.ny = ny;
    this.nr = nr;
    this.np = np;
    this.nk = nk;
    
    this.sym.rfun = this.sym.root;
    
    % initial hashes
    
    this.HTable.y = DataHash(char(this.sym.y));
    this.HTable.x = DataHash(char(this.sym.x));
    this.HTable.p = DataHash(char(this.sym.p));
    this.HTable.k = DataHash(char(this.sym.k));
    this.HTable.x0 = DataHash(char(this.sym.x0));
    this.HTable.rfun = DataHash(char(this.sym.rfun));
    if(strcmp(this.wtype,'iw'))
        this.HTable.xdot = DataHash(char(this.sym.xdot));
        this.HTable.dx0 = DataHash(char(this.sym.dx0));
        this.HTable.M = DataHash(char(this.sym.M));
    else
        this.HTable.xdot = DataHash(char(this.sym.xdot));
    end
    this.HTable.sigma_y = DataHash(char(this.sym.sigma_y));
    this.HTable.sigma_t = DataHash(char(this.sym.sigma_t));
    
    % compute functions
    
    % do not change the ordering, it is essential for correct dependencies
    funs = {'xdot','J','x0','Jv','JBand','JSparse','y','dydp','root','rootval','deltadisc','dxdotdp'};
    
    if(this.forward)
        funs = {funs{:},'sxdot','sx0','sy','sroot','s2root','srootval','s2rootval','sdeltadisc'};
    end
    if(this.adjoint)
        funs = {funs{:},'xBdot','qBdot','JB','JvB','JBandB','JSparseB','dydx','dtdx','drvaldx','dtdp','drvaldp','bdeltadisc','ideltadisc','sigma_y','sigma_t','dsigma_ydp','dsigma_tdp','sx0'};
    end
    
    if(strcmp(this.wtype,'iw'))
        funs = {funs{:},'dx0','sdx0'};
    end
    
    funs = unique(funs);
    
    this.funs = funs;
    
    % compute symbolic expressions
    for ifun = 1:length(funs)
        this = this.getFun(HTable,funs{ifun});
    end
    
    if(this.fun.J)
        fprintf('sparse | ')
        M = double(logical(this.sym.J~=sym(zeros(size(this.sym.J)))));
        this.sparseidx = find(M);
        
        [ubw,lbw] = ami_bandwidth(M);
        
        this.ubw = ubw;
        this.lbw = lbw;
        this.nnz = length(find(M(:)));
        
        I = arrayfun(@(x) find(M(:,x))-1,1:nx,'UniformOutput',false);
        this.rowvals = [];
        for ix = 1:nx
            this.colptrs(ix) = length(this.rowvals);
            this.rowvals = [this.rowvals; I{ix}];
        end
        this.colptrs(ix+1) = length(this.rowvals);

        if(this.adjoint)
            if(this.fun.JB)
                fprintf('sparseB | ')
                MB = transpose(M);
                this.sparseidxB = find(MB);
                I = arrayfun(@(x) find(MB(:,x))-1,1:nx,'UniformOutput',false);
                this.rowvalsB = [];
                for ix = 1:nx
                    this.colptrsB(ix) = length(this.rowvalsB);
                    this.rowvalsB = [this.rowvalsB; I{ix}];
                end
                this.colptrsB(ix+1) = length(this.rowvalsB);
            end
        end
    end

    % save hashtable

    HTable = this.HTable;
    nxtrue = this.nxtrue;
    nytrue = this.nytrue;
    nx = this.nx;
    ny = this.ny;
    np = this.np;
    nk = this.nk;
    ndisc = this.ndisc;
    nr = this.nr;
    nnonzeros = this.nnz;
    id = this.id;
    ubw = this.ubw;
    lbw = this.lbw;
    colptrs = this.colptrs;
    rowvals = this.rowvals;
    sparseidx = this.sparseidx;
    colptrsB = this.colptrsB;
    rowvalsB = this.rowvalsB;
    sparseidxB = this.sparseidxB;

    save(fullfile(this.wrap_path,'models',this.modelname,'hashes.mat'),'HTable','nxtrue','nytrue','nx','ny','np','nk','ndisc','nr','nnonzeros','id','ubw','lbw','colptrs','rowvals','sparseidx','colptrsB','rowvalsB','sparseidxB');

    fprintf('\r')

end

function [ubw,lbw] = ami_bandwidth(M)
    % ami_bandwidth implements a bandwidth function for older versionsn of MATLAB
    %
    % Parameters:
    %  M: matrix for which the bandwidth is to be computed
    %
    % Return values:
    %  ubw: upper matrix bandwidth
    %  lbw: lower matrix bandwidth
    [i,j] = find(M);
    ubw = max(max(j-i),0);
    lbw = max(max(i-j),0);
end