function [this,model] = getSyms(this,model)
    % getSyms computes the symbolic expression for the requested function
    %
    % Parameters:
    %  model: model definition object @type amimodel
    %
    % Return values:
    %  this: updated function definition object @type amifun
    %  model: updated model definition object @type amimodel
    
    % store often used variables for ease of notation, dependencies should
    % ensure that these variables are defined
    
    persistent x p sx w ndw
    
    nx = model.nx;
    nevent = model.nevent;
    np = model.np;
    nk = model.nk;
    nz = model.nz;
    ny = model.ny;
    
    fprintf([this.funstr ' | '])
    switch(this.funstr)
        case 'x'
            % create cell array of same size
            xs = cell(nx,1);
            % fill cell array
            for j=1:nx
                xs{j} = sprintf('x_%i',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(xs);
            x = this.sym;
            
        case 'dx'
            % create cell array of same size
            dxs = cell(nx,1);
            % fill cell array
            for j=1:nx
                dxs{j} = sprintf('dx_%i',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(dxs);
            
        case 'p'
            % create cell array of same size
            ps = cell(np,1);
            % fill cell array
            for j=1:np
                ps{j} = sprintf('p_%i',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(ps);
            p = this.sym;
            
        case 'k'
            % create cell array of same size
            ks = cell(nk,1);
            % fill cell array
            for j=1:nk
                ks{j} = sprintf('k_%i',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(ks);
            
        case 'sx'
            % create cell array of same size
            sxs = cell(nx,np);
            % fill cell array
            for j = 1:nx
                for i = 1:np
                    sxs{j,i} = sprintf('sx_%i', j-1);
                end
            end
            % transform into symbolic expression
            this.sym = sym(sxs);
            sx = this.sym;
            
        case 'sdx'
            % create cell array of same size
            sdx = cell(nx,np);
            % fill cell array
            for j = 1:nx
                for i = 1:np
                    sdx{j,i} = sprintf('sdx_%i', j-1);
                end
            end
            % transform into symbolic expression
            this.sym = sym(sdx);
            
        case 'xB'
            % create cell array of same size
            xBs = cell(nx,1);
            % fill cell array
            for j = 1:nx
                xBs{j} = sprintf('xB_%i', j-1);
            end
            % transform into symbolic expression
            this.sym = sym(xBs);
            
        case 'dxB'
            % create cell array of same size
            dxBs = cell(nx,1);
            % fill cell array
            for j = 1:nx
                dxBs{j} = sprintf('dxB_%i', j-1);
            end
            % transform into symbolic expression
            this.sym = sym(dxBs);
            
        case 'y'
            this.sym = model.sym.y;
            % replace unify symbolic expression
            this = unifySyms(this,model);
            % create cell array of same size
            ys = cell(ny,1);
            % fill cell array
            for j = 1:ny
                ys{j} = sprintf('y_%i', j-1);
            end
            % transform into symbolic expression
            this.strsym = sym(ys);
            
            
            % activate splines
            for iy = 1:ny
                if(not(all([model.splineflag,model.minflag,model.maxflag])))
                    str = char(this.sym(iy));
                    if(strfind(str,'spline'))
                        model.splineflag = true;
                    end
                    if(strfind(str,'max'))
                        model.maxflag = true;
                    end
                    if(strfind(str,'min'))
                        model.minflag = true;
                    end
                end
            end
            
            
        case 'x0'
            this.sym = model.sym.x0;
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'dx0'
            this.sym = model.sym.dx0;
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'sigma_y'
            this.sym = model.sym.sigma_y;
            sdys = cell(model.nytrue,1);
            % fill cell array
            for j=1:model.nytrue
                sdys{j} = sprintf('sdy_%i',j-1);
            end
            % transform into symbolic expression
            this.strsym = sym(sdys);
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'sigma_z'
            this.sym = model.sym.sigma_z;
            sdzs = cell(model.nztrue,1);
            % fill cell array
            for j=1:model.nztrue
                sdzs{j} = sprintf('sdz_%i',j-1);
            end
            % transform into symbolic expression
            this.strsym = sym(sdzs);
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'M'
            this.sym = model.sym.M;
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'xdot'
            this.sym = simplify(model.sym.xdot);
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
            if(strcmp(model.wtype,'iw'))
                if(size(this.sym,2)>size(this.sym,1))
                    this.sym = -transpose(model.fun.M.sym*model.fun.dx.sym)+this.sym;
                else
                    this.sym = -model.fun.M.sym*model.fun.dx.sym+this.sym;
                end
            end
            
            % create cell array of same size
            xdots = cell(nx,1);
            xdot_olds = cell(nx,1);
            % fill cell array
            for j=1:nx
                xdots{j} = sprintf('xdot_%i',j-1);
                xdot_olds{j} = sprintf('xdot_old_%i',j-1);
            end
            this.strsym = sym(xdots);
            this.strsym_old = sym(xdot_olds);
            
            % activate splines
            for ix = 1:nx
                if(not(all([model.splineflag,model.minflag,model.maxflag])))
                    str = char(this.sym(ix));
                    if(strfind(str,'spline'))
                        model.splineflag = true;
                    end
                    if(strfind(str,'max'))
                        model.maxflag = true;
                    end
                    if(strfind(str,'min'))
                        model.minflag = true;
                    end
                end
            end
            
        case 'w'
            optimize = getoptimized(optsym(model.fun.xdot.sym));
            tmpxdot = children(optimize(end));
            temps = sym('t',[(length(optimize)),1]);
            nw = (length(optimize)-1);
            model.nw = nw;
            exprs = arrayfun(@(x) children(x),optimize(1:(end-1)),'UniformOutput',false); % extract symbolic variable
            S.subs = {2};
            S.type='()';
            C = cellfun(@(x) subsref(x,S),exprs,'UniformOutput',false); % get second element
            this.sym = [C{:}]; % transform cell to matrix
%             model.nw = 0;
%             nw = 0;
%             this.sym = sym(zeros(0,1));          
            

            
            ws = cell(nw,1);
            % fill cell array
            for iw = 1:nw
                ws{iw} = sprintf('w_%i', iw-1);
            end
            % transform into symbolic expression
            this.strsym = sym(ws);
            tmpxdot = mysubs(tmpxdot,temps(2:end),this.strsym); % replace common expressions
            this.sym = mysubs(this.sym,temps(2:end),this.strsym);
            model.updateRHS(tmpxdot); % update rhs
            
            w = this.strsym;
            
            % find hierarchy depth
            ndw = 0;
            jacw = jacobian(model.fun.w.sym,w);
            S = sparse(double(jacw~=0));
            Sndw = S;
            while(any(any(Sndw)))
                ndw = ndw+1;
                Sndw = Sndw*S;
            end
         
        case 'dwdx'
            jacx = jacobian(model.fun.w.sym,x);
            jacw = jacobian(model.fun.w.sym,w);
            this.sym = jacx + jacw*jacx;
            for ind = 2:ndw
                this.sym = this.sym + jacw^ind*jacx;
            end
            % fill cell array
            idx_w = find(this.sym~=0);
            dwdxs = cell(model.nw,nx);
            [dwdxs{:}] = deal('0');
            if(numel(idx_w)>0)
                for iw = transpose(idx_w)
                    dwdxs{iw} = sprintf('dwdx_%i', iw-1);
                end
            end
            % transform into symbolic expression
            this.strsym = sym(dwdxs);
            
        case 'dwdp'
            jacp = jacobian(model.fun.w.sym,p);
            jacw = jacobian(model.fun.w.sym,w);
            this.sym = jacp + jacw*jacp;
            for ind = 2:ndw
                this.sym = this.sym + jacw^ind*jacp;
            end
            this.sym = jacobian(model.fun.w.sym,p);
            % fill cell array
            idx_w = find(this.sym~=0);
            dwdps = cell(model.nw,np);
            [dwdps{:}] = deal('0');
            if numel(idx_w)>0
                for iw = transpose(idx_w)
                    dwdps{iw} = sprintf('dwdx_%i', iw-1);
                end
            end
            % transform into symbolic expression
            this.strsym = sym(dwdps);
            
        case 'dfdx'
            this.sym=jacobian(model.fun.xdot.sym,x) + jacobian(model.fun.xdot.sym,w)*model.fun.dwdx.strsym;
            
        case 'J'
            if(strcmp(model.wtype,'iw'))
                syms cj
                this.sym = model.fun.dfdx.sym - cj*model.fun.M.sym;
            else
                this.sym = jacobian(model.fun.xdot.sym,x)  + jacobian(model.fun.xdot.sym,w)*model.fun.dwdx.strsym;
            end
            
            %%
            % build short strings for reuse of jacobian
            
            % find nonzero entries
            idx = find(logical(this.sym~=0));
            % create cell array of same size
            Js = cell(length(idx),1);
            % fill cells with strings
            for iJ = 1:length(idx)
                Js{iJ} = sprintf('tmp_J%i',iJ-1);
            end
            % create full symbolic matrix
            this.strsym = sym(zeros(nx,nx));
            % fill nonzero entries
            this.strsym(idx) = sym(Js);
            
        case 'JB'
            this.sym=-transpose(model.fun.J.sym);
            
        case 'dxdotdp'
            this.sym=jacobian(model.fun.xdot.sym,p)  + jacobian(model.fun.xdot.sym,w)*model.fun.dwdp.strsym;
            
            %%
            % build short strings for reuse of dxdotdp
            
            % find nonzero rows, only completely nonzero rows matter as we
            % will have parameter dependent implementation later on
            % anyways.
            idx = find(any(logical(this.sym~=0),2));
            % create cell array of same size
            dxdotdps = cell(length(idx),1);
            % fill cells with strings
            for ix = 1:length(idx)
                dxdotdps{ix} = sprintf('tmp_dxdotdp%i',ix-1);
            end
            % create full symbolic array
            this.strsym = sym(zeros(nx,1));
            % fill nonzero entries
            this.strsym(idx) = sym(dxdotdps);
            
        case 'sx0'
            this.sym=jacobian(model.fun.x0.sym,p);
            
        case 'sdx0'
            this.sym=jacobian(model.fun.dx0.sym,p);
            
        case 'sxdot'
            if(strcmp(model.wtype,'iw'))
                this.sym=model.fun.dfdx.sym*sx-model.fun.M.sym*model.fun.sdx.sym+model.fun.dxdotdp.sym;
            else
                this.sym=model.fun.J.strsym*sx(:,1)+model.fun.dxdotdp.strsym;
            end
            
        case 'dydx'
            this.sym=jacobian(model.fun.y.sym,x);
            % create cell array of same size
            dydxs = cell(ny,nx);
            % fill cell array
            for j = 1:ny
                for i = 1:nx
                    if(this.sym(j,i)~=0)
                        dydxs{j,i} = sprintf('dydx_%i', j-1 + (i-1)*ny);
                    else
                        dydxs{j,i} = sprintf('0');
                    end
                end
            end
            % transform into symbolic expression
            this.strsym = sym(dydxs);
            
        case 'dydp'
            this.sym=jacobian(model.fun.y.sym,p);
            % create cell array of same size
            dydps = cell(ny,np);
            % fill cell array
            for j = 1:ny
                for i = 1:np
                    if(this.sym(j,i)~=0)
                        dydps{j,i} = sprintf('dydp_%i', j-1 + (i-1)*ny);
                    else
                        dydps{j,i} = sprintf('0');
                    end
                end
            end
            % transform into symbolic expression
            this.strsym = sym(dydps);
            
        case 'sy'
            this.sym=model.fun.dydp.sym + model.fun.dydx.sym*model.fun.sx.sym ;
            % create cell array of same size
            sys = cell(ny,np);
            % fill cell array
            for j = 1:ny
                for i = 1:np
                    sys{j,i} = sprintf('sy_%i', j-1);
                end
            end
            % transform into symbolic expression
            this.strsym = sym(sys);
            
        case 'Jv'
            % create cell array of same size
            vs = cell(nx,1);
            % fill cells
            for j=1:nx
                vs{j} = sprintf('v_%i',j-1);
            end
            % transform to symbolic variable
            vs = sym(vs);
            % multiply
            this.sym = model.fun.J.sym*vs;
            
        case 'JvB'
            % create cell array of same size
            vs = cell(nx,1);
            % fill cells
            for j=1:nx
                vs{j} = sprintf('v_%i',j-1);
            end
            % transform to symbolic variable
            vs = sym(vs);
            % multiply
            this.sym = -transpose(model.fun.J.sym)*vs;
            
        case 'xBdot'
            if(strcmp(model.wtype,'iw'))
                syms t
                this.sym = diff(transpose(model.fun.M.sym),t)*model.fun.xB.sym + transpose(model.fun.M.sym)*model.fun.dxB.sym - transpose(model.fun.dfdx.sym)*model.fun.xB.sym;
            else
                this.sym = -transpose(model.fun.J.sym)*model.fun.xB.sym;
            end
            
        case 'qBdot'
            this.sym = -transpose(model.fun.xB.sym)*model.fun.dxdotdp.strsym;
            
        case 'dsigma_ydp'
            this.sym = jacobian(model.fun.sigma_y.sym,p);
            % create cell array of same size
            dsdydps = cell(model.nytrue,np);
            % fill cell array
            for iy = 1:model.nytrue
                for ip = 1:np
                    if(this.sym(iy,ip)~=0)
                        dsdydps{iy,ip} = sprintf('dsdydp_%i', iy-1 + (ip-1)*ny);
                    else
                        dsdydps{iy,ip} = sprintf('0');
                    end
                end
            end
            % transform into symbolic expression
            this.strsym = sym(dsdydps);
            
        case 'dsigma_zdp'
            if(nz>0)
                this.sym = jacobian(model.fun.sigma_z.sym,p);
            else
                this.sym = sym(zeros(model.nztrue,np));
            end
            % create cell array of same size
            dsdydps = cell(model.nztrue,np);
            % fill cell array
            for j = 1:model.nztrue
                for i = 1:np
                    if(this.sym(j,i)~=0)
                        dsdydps{j,i} = sprintf('dsdzdp_%i', j-1 + (i-1)*nz);
                    else
                        dsdydps{j,i} = sprintf('0');
                    end
                end
            end
            % transform into symbolic expression
            this.strsym = sym(dsdydps);
            
        case 'root'
            if(nevent>0)
                this.sym = transpose([model.event.trigger]);
                this = unifySyms(this,model);
            else
                this.sym = sym(zeros(0,1));
            end
            
        case 'drootdp'
            this.sym = jacobian(model.fun.root.sym,p);
            
        case 'drootdx'
            this.sym = jacobian(model.fun.root.sym,x);
            
        case 'drootdt'
            this.sym = diff(model.fun.root.sym,'t') + model.fun.drootdx.sym*model.fun.xdot.sym;
            
        case 'sroot'
            this.sym = model.fun.drootdp.sym + model.fun.drootdx.sym*sx;
            
        case 'dtaudp'
            this.sym = sym(zeros(nevent,np));
            for ievent = 1:nevent
                this.sym(ievent,:) = - model.fun.drootdp.sym(ievent,:)/model.fun.drootdt.sym(ievent);
            end
            
        case 'stau'
            this.sym = sym(zeros(nevent,np));
            for ievent = 1:nevent
                this.sym(ievent,:) = - model.fun.sroot.sym(ievent,:)/model.fun.drootdt.sym(ievent);
            end
            
        case 'deltax'
            if(nevent>0)
                this.sym = [model.event.bolus];
                this = unifySyms(this,model);
            else
                this.sym = sym(zeros(0,1));
            end
            
        case 'deltaxdot'
            this.sym = model.fun.xdot.strsym-model.fun.xdot.strsym_old;
            
        case 'ddeltaxdp'
            this.sym = sym(zeros(nx,nevent,np));
            for ievent = 1:nevent
                this.sym(:,ievent,:) = jacobian(model.fun.deltax.sym(:,ievent),p);
            end
            
        case 'ddeltaxdx'
            this.sym = sym(zeros(nx,nevent,nx));
            for ievent = 1:nevent
                this.sym(:,ievent,:) = jacobian(model.fun.deltax.sym(:,ievent),x);
            end
            
        case 'ddeltaxdt'
            this.sym = diff(model.fun.deltax.sym,'t');
            
        case 'deltasx'
            
            if(nevent>0)
                for ievent = 1:nevent
                    
                    % dtdp  = (1/drdt)*drdp
                    dtdp = model.fun.stau.sym(ievent,:);
                    
                    % if we are just non-differentiable and but not
                    % discontinuous we can ignore some of the terms!                
                    if(any(logical(model.fun.deltax.sym(:,ievent)~=0)))
                        % dxdp  = dx/dt*dt/dp + dx/dp
                        dxdp = model.fun.xdot.sym*dtdp + sx;
                        
                        this.sym(:,:,ievent) = ...
                            + permute(model.fun.ddeltaxdx.sym(:,ievent,:),[1 3 2])*dxdp ...
                            + model.fun.ddeltaxdt.sym(:,ievent)*dtdp ...
                            + permute(model.fun.ddeltaxdp.sym(:,ievent,:),[1 3 2]);
                    else
                        this.sym(:,:,ievent) = sym(zeros([nx,np]));
                    end
                    if(any(model.event(ievent).hflag))
                        this.sym(model.event(ievent).hflag,:,ievent) = ...
                            this.sym(model.event(ievent).hflag,:,ievent) ...
                            - model.fun.deltaxdot.sym(model.event(ievent).hflag)*dtdp;
                    end
                end
            end
            
        case 'deltaqB'
            this.sym = sym(zeros(np,nevent));
            for ievent = 1:nevent
                this.sym(:,ievent) = transpose(model.fun.xB.sym)*squeeze(model.fun.ddeltaxdp.sym(:,ievent,:));
            end
            
        case 'deltaxB'
            this.sym = sym(zeros(nx,nevent));
            for ievent = 1:nevent
                this.sym(:,ievent) = squeeze(model.fun.ddeltaxdx.sym(:,ievent,:))*model.fun.xB.sym;
            end
            
        case 'z'
            if(nevent>0)
                this.sym = transpose([model.event.z]);
                this = unifySyms(this,model);
            else
                this.sym = sym(zeros(0,1));
            end
            % construct the event identifyer, this is a vector which maps
            % events to outputs z
            model.z2event = zeros(nz,1);
            iz = 0;
            for ievent = 1:nevent
                for jz = 1:length(model.event(ievent).z)
                    iz = iz+1;
                    model.z2event(iz) = ievent;
                end
            end
            % create cell array of same size
            zs = cell(nz,1);
            % fill cell array
            for j = 1:nz
                zs{j} = sprintf('z_%i', j-1);
            end
            % transform into symbolic expression
            this.strsym = sym(zs);
            
        case 'dzdp'
            this.sym = jacobian(model.fun.z.sym,p);
            % create cell array of same size
            dzdps = cell(nz,np);
            % fill cell array
            for j = 1:nz
                for i = 1:np
                    if(this.sym(j,i)~=0)
                        dzdps{j,i} = sprintf('dzdp_%i', j-1 + (i-1)*nz);
                    else
                        dzdps{j,i} = sprintf('0');
                    end
                end
            end
            % transform into symbolic expression
            this.strsym = sym(dzdps);
            
        case 'dzdx'
            this.sym = jacobian(model.fun.z.sym,x);
            % create cell array of same size
            dzdxs = cell(nz,nx);
            % fill cell array
            for j = 1:nz
                for i = 1:nx
                    if(this.sym(j,i)~=0)
                        dzdxs{j,i} = sprintf('dzdx_%i', j-1 + (i-1)*nz);
                    else
                        dzdxs{j,i} = sprintf('0');
                    end
                end
            end
            % transform into symbolic expression
            this.strsym = sym(dzdxs);
            
        case 'dzdt'
            this.sym = diag(diff(model.fun.z.sym,'t'));
            
        case 'sz'
            % the tau is event specific so mapping from z to events is
            % necessary here.
            this.sym = ...
                + model.fun.dzdp.sym ...
                + model.fun.dzdx.sym*sx ...
                + model.fun.dzdt.sym*model.fun.stau.sym(model.z2event,:);
            % create cell array of same size
            szs = cell(nz,np);
            % fill cell array
            for j = 1:nz
                for i = 1:np
                    szs{j,i} = sprintf('sz_%i', j-1);
                end
            end
            % transform into symbolic expression
            this.strsym = sym(szs);
            
        case 'sz_tf'
            % the tau is event specific so mapping from z to events is
            % necessary here.
            this.sym = ...
                + model.fun.dzdp.sym ...
                + model.fun.dzdx.sym*sx;
            
        case 'JBand'
            %do nothing
        case 'JBandB'
            %do nothing
        case 'JSparse'
            %do nothing
        case 'JSparseB'
            %do nothing
            
        case 'Jy'
            this.sym = model.sym.Jy(:);
        case 'dJydy'
            this.sym = jacobian(model.fun.Jy.sym,model.fun.y.strsym);
        case 'dJydx'
            this.sym = model.fun.dJydy.sym*model.fun.dydx.strsym;
        case 'dJydsigma'
            this.sym = jacobian(model.fun.Jy.sym,model.fun.sigma_y.strsym);
        case 'dJydp'
            this.sym = model.fun.dJydy.sym*model.fun.dydp.strsym + model.fun.dJydsigma.sym*model.fun.dsigma_ydp.strsym;
        case 'sJy'
            this.sym = model.fun.dJydy.sym*model.fun.sy.strsym + model.fun.dJydp.sym;
            
        case 'Jz'
            this.sym = model.sym.Jz(:);
        case 'dJzdz'
            if(nz>0)
                this.sym = jacobian(model.fun.Jz.sym,model.fun.z.strsym);
            else
                this.sym = sym(zeros(1,0));
            end
        case 'dJzdx'
            this.sym = model.fun.dJzdz.sym*model.fun.dzdx.strsym;
        case 'dJzdsigma'
            if(nz>0)
                this.sym = jacobian(model.fun.Jz.sym,model.fun.sigma_z.strsym);
            else
                this.sym = sym(zeros(1,0));
            end
        case 'dJzdp'
            this.sym = model.fun.dJzdz.sym*model.fun.dzdp.strsym + model.fun.dJzdsigma.sym*model.fun.dsigma_zdp.strsym;
        case 'sJz'
            this.sym = model.fun.dJzdz.sym*model.fun.sz.strsym + model.fun.dJzdp.sym;
            
            
        otherwise
            error('unknown function name')
    end
end

function this = unifySyms(this,model)
    % unifySyms replaces model specific variable names with general
    % ones which conforms with C naming schemes
    %
    % Parameters:
    %  model: model definition object @type amimodel
    this.sym = mysubs(this.sym,model.sym.x,model.fun.x.sym);
    this.sym = mysubs(this.sym,model.sym.p,model.fun.p.sym);
    this.sym = mysubs(this.sym,model.sym.k,model.fun.k.sym);
end

function out = mysubs(in, old, new)
    % mysubs is a wrapper for ther subs matlab function
    %
    % Parameters:
    %  in: symbolic expression in which to replace @type symbolic
    %  old: expression to be replaced @type symbolic
    %  new: replacement expression @type symbolic
    %
    % Return values:
    %  out: symbolic expression with replacement @type symbolic
    if(~isnumeric(in) && ~isempty(old) && ~isempty(symvar(in)))
        matVer = ver('MATLAB');
        if(str2double(matVer.Version)>=8.1)
            out = subs(in, old(:), new(:));
        else
            out = subs(in, old(:), new(:), 0);
        end
    else
        out = in;
    end
end