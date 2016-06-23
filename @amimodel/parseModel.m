function parseModel(this)
    % parseModel parses the model definition and computes all necessary symbolic expressions.
    %
    
    % compile CalcMD5 if necessary
    try
        CalcMD5('TEST','char','hex');
    catch
        try
            addpath(fullfile(this.wrap_path,'auxiliary','CalcMD5'))
            CalcMD5('TEST','char','hex');
        catch
            disp('CalcMD5 has not been compiled yet. Compiling now!')
            tmpdir = pwd;
            cd(fullfile(this.wrap_path,'auxiliary','CalcMD5'))
            mex(fullfile(this.wrap_path,'auxiliary','CalcMD5','CalcMD5.c'))
            addpath(fullfile(this.wrap_path,'auxiliary','CalcMD5'))
            cd(tmpdir);
        end
    end

    % load old hashes
    HTable = this.loadOldHashes();
    
    nevent = length(this.event);
    this.nevent = nevent;
    
    % fix z2event if z was computed previously
    if(isfield(this.fun,'z'))
        this.fun = rmfield(this.fun,'z');
        this.getFun([],'z');
    end
    
    nx = length(this.sym.x);
    np = length(this.sym.p);
    nk = length(this.sym.k);
    ny = length(this.sym.y);
    nz = length([this.event.z]);
    if(this.nytrue == 0) % only set this if it still has the default value, if 0 is already the non default value it does not matter anyways
        nytrue = length(this.sym.sigma_y);
        this.nytrue = nytrue;
    end  
    if(this.nztrue == 0)
        nztrue = length(this.sym.sigma_z);
        this.nztrue = nztrue;
    end
    
    %check zero trigger events
    for ievent = 1:nevent
        if(isequaln(this.event(ievent).trigger,sym(0)))
            error('Trigger functions cannot be equal to zero. Please check your event definition!')
        end
    end
    
    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.np = np;
    this.nk = nk;
    
    % initial hashes
    this.HTable(1).y = CalcMD5(char(this.sym.y));
    this.HTable(1).x = CalcMD5(char(this.sym.x));
    this.HTable(1).p = CalcMD5(char(this.sym.p));
    this.HTable(1).k = CalcMD5(char(this.sym.k));
    this.HTable(1).x0 = CalcMD5(char(this.sym.x0));
    if(nevent>0)
        this.HTable(1).trigger = CalcMD5(char([this.event.trigger]));
        this.HTable(1).bolus = CalcMD5(char([this.event.bolus]));
        this.HTable(1).z = CalcMD5(char([this.event.z]));
    end
    if(strcmp(this.wtype,'iw'))
        this.HTable(1).xdot = CalcMD5(char(this.sym.xdot));
        this.HTable(1).dx0 = CalcMD5(char(this.sym.dx0));
        this.HTable(1).M = CalcMD5(char(this.sym.M));
    else
        this.HTable(1).xdot = CalcMD5(char(this.sym.xdot));
    end
    this.HTable(1).sigma_y = CalcMD5(char(this.sym.sigma_y));
    this.HTable(1).sigma_z = CalcMD5(char(this.sym.sigma_z));
    
    % compute functions
    
    % do not change the ordering, it is essential for correct dependencies
    funs = {'xdot','w','dwdx','dwdp','J','x0','Jv','JBand','JSparse','y','z','deltax','dydp','dxdotdp','root','Jy','dJydx','dJydp','sJy','Jz','dJzdx','dJzdp','sJz'};
    
    if(this.forward)
        funs = {funs{:},'sxdot','sx0','sy','sz','sz_tf','deltasx','stau','sroot'};
    end
    if(this.adjoint)
        funs = {funs{:},'xBdot','qBdot','JB','JvB','JBandB','JSparseB','dydx','dzdx','dzdp','deltaxB','deltaqB','sigma_y','sigma_z','dsigma_ydp','dsigma_zdp','sx0'};
    end
    if(this.o2flag)
        funs = {funs{:},'s2root'};
    end
    
    if(strcmp(this.wtype,'iw'))
        funs = {funs{:},'dx0','sdx0','dfdx','M'};
    end
    
    funs = unique(funs);
    
    this.funs = funs;
    
    % compute symbolic expressions
    for ifun = 1:length(funs)
        this.getFun(HTable,funs{ifun});
    end
    
    if(isfield(this.fun,'J'))
        fprintf('sparse | ')
        M = double(logical(this.fun.J.sym~=sym(zeros(size(this.fun.J.sym)))));
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
            if(isfield(this.fun,'JB'))
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
    
    if(strcmp(this.wtype,'iw'))
        this.getFun([],'M');
        this.id = double(logical(sum(this.fun.M.sym,2)~=0));
    else
        this.id = zeros(nx,1);
    end

    % save hashtable
    HTable = this.HTable;
    nxtrue = this.nxtrue;
    nytrue = this.nytrue;
    nx = this.nx;
    ny = this.ny;
    np = this.np;
    nk = this.nk;
    nz = this.nz;
    nw = this.nw;
    ndwdx = this.ndwdx;
    ndwdp = this.ndwdp;
    nevent = this.nevent;
    z2event = this.z2event;
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
    compver = this.compver;

    save(fullfile(this.wrap_path,'models',this.modelname,'hashes.mat'),'HTable','nxtrue','nytrue','nx','ny','np','nk','nevent','nz','z2event','nnonzeros','id','ubw','lbw','colptrs','rowvals','sparseidx','colptrsB','rowvalsB','sparseidxB','compver','ndwdx','ndwdp','nw');

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
    if(isempty(M))
        ubw = 0;
        lbw = 0;
    else
        [i,j] = find(M);
        ubw = max(max(j-i),0);
        lbw = max(max(i-j),0);
    end
end