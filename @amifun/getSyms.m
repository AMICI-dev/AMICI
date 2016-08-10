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
    
    persistent x p sx w ndw jacw
    
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
            sxs = cell(nx,1);
            % fill cell array
            for j = 1:nx
                sxs{j} = sprintf('sx_%i', j-1);
            end
            % transform into symbolic expression
            this.sym = repmat(sym(sxs),[1,np]);
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
            this = makeStrSymsFull(this);
            
            
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
            this = makeStrSyms(this);
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'sigma_z'
            this.sym = model.sym.sigma_z;
            this = makeStrSyms(this);
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'M'
            this.sym = model.sym.M;
            % replace unify symbolic expression
            this = unifySyms(this,model);
            this = makeStrSyms(this);
            
        case 'xdot'
            this.sym = model.sym.xdot;
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
            tmpxdot = sym(char(optimize(end)));
            nw = (length(optimize)-1);
            model.nw = nw;
            if(nw>0)
                exprs = arrayfun(@(x) children(x),optimize(1:(end-1)),'UniformOutput',false); % extract symbolic variable
                S.subs = {2};
                S.type='()';
                C = cellfun(@(x) subsref(x,S),exprs,'UniformOutput',false); % get second element
                this.sym = [C{:}]; % transform cell to matrix
                S.subs = {1};
                C = cellfun(@(x) subsref(x,S),exprs,'UniformOutput',false);
                temps = [C{:}];
            end
%             model.nw = 0;
%             nw = 0;
%             this.sym = sym(zeros(0,1));          
            

            
            ws = cell(nw,1);
            ts = cell(nw,1);
            % fill cell array
            for iw = 1:nw
                ws{iw} = sprintf('w_%i', iw-1);
            end
            % transform into symbolic expression
            this.strsym = sym(ws);
            ndw = 0;
            if(nw>0)
                tmpxdot = mysubs(tmpxdot,temps,this.strsym); % replace common expressions
                this.sym = mysubs(this.sym,temps,this.strsym);
                model.updateRHS(tmpxdot); % update rhs
                ndw = 1;
                jacw = jacobian(this.sym,this.strsym);
                vv = sym('v',[length(this.strsym),1]);
                while(sum(jacw^ndw*vv)~=0)
                    ndw = ndw+1;
                end
                ndw = ndw - 1;
            end
            w = this.strsym;

        case 'dwdx'
            if(length(model.fun.w.sym)>0)
                jacx = jacobian(model.fun.w.sym,x);
                this.sym = jacx;
                for idw = 1:ndw
                    this.sym = this.sym + (jacw^idw)*jacx; % this part is only to get the right nonzero entries 
                end
                % fill cell array
                idx_w = find(this.sym);
                this.strsym = sym(zeros(size(jacx)));
                if(numel(idx_w)>0)
                    for iw = 1:length(idx_w)
                        this.strsym(idx_w(iw)) = sym(sprintf('dwdx_%i', iw-1));
                    end
                    model.ndwdx = length(idx_w);
                    % update dwdx with simplified expressions, here we can exploit
                    % the proper ordering of w to ensure correctness of expressions
                    tmp = jacx + jacw*this.strsym;
                    this.sym = tmp(idx_w);
                else
                    this.strsym = sym(zeros(size(jacx)));
                end
            else
                this.sym = sym(zeros(0,nx));
                this.strsym = sym(zeros(0,nx));
            end
            
        case 'dwdp'
            if(length(model.fun.w.sym)>0)
                jacp = jacobian(model.fun.w.sym,p);
                this.sym = jacp;
                for idw = 1:ndw
                    this.sym = this.sym + (jacw^idw)*jacp; % this part is only to get the right nonzero entries 
                end
                % fill cell array
                idx_w = find(this.sym);
                this.strsym = sym(zeros(size(jacp)));
                if(numel(idx_w)>0)
                    for iw = 1:length(idx_w)
                        this.strsym(idx_w(iw)) = sym(sprintf('dwdp_%i', iw-1));
                    end
                    model.ndwdp = length(idx_w);
                    % update dwdx with simplified expressions, here we can exploit
                    % the proper ordering of w to ensure correctness of expressions
                    tmp = jacp + jacw*this.strsym;
                    this.sym = tmp(idx_w);
                end
            else
                this.sym = sym(zeros(0,nx));
                this.strsym = sym(zeros(0,nx));
            end
            
        case 'dfdx'
            this.sym=jacobian(model.fun.xdot.sym,x) + jacobian(model.fun.xdot.sym,w)*model.fun.dwdx.strsym;
            this = makeStrSyms(this);
            
        case 'J'
            if(strcmp(model.wtype,'iw'))
                syms cj
                this.sym = model.fun.dfdx.sym - cj*model.fun.M.sym;
            else
                if(~isempty(w))
                    this.sym = jacobian(model.fun.xdot.sym,x)  + jacobian(model.fun.xdot.sym,w)*model.fun.dwdx.strsym;
                    this.sym_noopt = jacobian(model.fun.xdot.sym_noopt,x);
                else
                    this.sym = jacobian(model.fun.xdot.sym,x);
                    this.sym_noopt = this.sym;
                end
            end
            
            %%
            % build short strings for reuse of jacobian
            
            % find nonzero entries
            idx = find(this.sym);
            % create cell array of same size
            Js = sym(zeros(length(idx),1));
            % fill cells with strings
            for iJ = 1:length(idx)
                Js(iJ) = sym(sprintf('tmp_J%i',iJ-1));
            end
            % create full symbolic matrix
            this.strsym = sym(zeros(nx,nx));
            % fill nonzero entries
            this.strsym(idx) = Js;
            
        case 'JB'
            this.sym = sym(zeros(model.nx, model.nx));
            % Augmentation needs a transposition in the submatrices
            for ix = 1 : model.ng
                for jx = 1 : model.ng
                    this.sym((ix-1)*model.nxtrue+1:ix*model.nxtrue, (jx-1)*model.nxtrue+1:jx*model.nxtrue) = ...
                        -transpose(model.fun.J.sym((ix-1)*model.nxtrue+1:ix*model.nxtrue, (jx-1)*model.nxtrue+1:jx*model.nxtrue));
                end
            end
            
        case 'dxdotdp'
            if(~isempty(w))
                this.sym=jacobian(model.fun.xdot.sym,p)  + jacobian(model.fun.xdot.sym,w)*model.fun.dwdp.strsym;
                this.sym_noopt = jacobian(model.fun.xdot.sym_noopt,p);
            else
                this.sym=jacobian(model.fun.xdot.sym,p);
                this.sym_noopt = this.sym;
            end
            
            %%
            % build short strings for reuse of dxdotdp
            % create cell array of same size
            dxdotdps = cell(nx,1);
            % fill cells with strings
            for ix = 1:nx
                dxdotdps{ix} = sprintf('tmp_dxdotdp%i',ix-1);
            end
            % create full symbolic array
            this.strsym = sym(dxdotdps);
            
        case 'sx0'
            this.sym=jacobian(model.fun.x0.sym,p);
            
        case 'sdx0'
            this.sym=jacobian(model.fun.dx0.sym,p);
            
        case 'sxdot'
            if(np>0)
                if(strcmp(model.wtype,'iw'))
                    this.sym=model.fun.dfdx.strsym*sx(:,1)-model.fun.M.strsym*model.fun.sdx.sym(:,1)+model.fun.dxdotdp.strsym;
                else
                    this.sym=model.fun.J.strsym*sx(:,1)+model.fun.dxdotdp.strsym;
                end
            else
                this.sym = sym(zeros(size(sx,1),0));
            end
            
        case 'dydx'
            this.sym=jacobian(model.fun.y.sym,x);
            % create cell array of same sizex
            this.strsym = sym(zeros(ny,nx));
            % fill cell array
            this = makeStrSyms(this);
            
        case 'dydp'
            this.sym=jacobian(model.fun.y.sym,p);
            % create cell array of same size
            this = makeStrSyms(this);
            
        case 'sy'
            this.sym=model.fun.dydp.strsym + model.fun.dydx.strsym*model.fun.sx.sym ;
            % create cell array of same size
            this = makeStrSyms(this);
            
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
            this.sym = -model.fun.JB.sym*vs;
            
        case 'xBdot'
            if(strcmp(model.wtype,'iw'))
                syms t
                this.sym = diff(transpose(model.fun.M.sym),t)*model.fun.xB.sym + transpose(model.fun.M.sym)*model.fun.dxB.sym - transpose(model.fun.dfdx.sym)*model.fun.xB.sym;
            else
                % Augmenting the system needs transposition of submatrices
                % I'm sure, there is a more intelligent solution to it...
                this.sym = sym(zeros(nx,1));
                if(model.o2flag)
                    for ix = 1 : model.ng
                        for jx = 1 : model.ng
                            this.sym((ix-1)*model.nxtrue+1 : ix*model.nxtrue) = this.sym((ix-1)*model.nxtrue+1 : ix*model.nxtrue) - ...
                                transpose(model.fun.J.sym((ix-1)*model.nxtrue+1 : ix*model.nxtrue , (jx-1)*model.nxtrue+1 : jx*model.nxtrue)) ...
                                * model.fun.xB.sym((jx-1)*model.nxtrue+1 : jx*model.nxtrue);
                        end
                    end
                else
                    this.sym = transpose(model.fun.J.sym) * model.fun.xB.sym;
                end
            end
            
        case 'qBdot'
            % If we do second order adjoints, we have to augment
            if (model.nxtrue < nx)
                this.sym = sym(zeros(model.ng, model.np));
                this.sym(1,:) = -transpose(model.fun.xB.sym(1:model.nxtrue)) * model.fun.dxdotdp.sym(1:model.nxtrue, :);
                for ig = 2 : model.ng
                    this.sym(ig,:) = ...
                        -transpose(model.fun.xB.sym(1:model.nxtrue)) * model.fun.dxdotdp.sym((ig-1)*model.nxtrue+1 : ig*model.nxtrue, :) ...
                        -transpose(model.fun.xB.sym((ig-1)*model.nxtrue+1 : ig*model.nxtrue)) * model.fun.dxdotdp.sym(1:model.nxtrue, :);
                end                    
            else
                this.sym = -transpose(model.fun.xB.sym)*model.fun.dxdotdp.sym;
            end

        case 'dsigma_ydp'
            this.sym = jacobian(model.fun.sigma_y.sym,p);
            this = makeStrSyms(this);
            
        case 'dsigma_zdp'
            if(nz>0)
                this.sym = jacobian(model.fun.sigma_z.sym,p);
            else
                this.sym = sym(zeros(model.nz,np));
            end
            this = makeStrSyms(this);
            
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
          
        case 's2root'
            switch(model.o2flag)
                case 1
                    s2x = reshape(sx((model.nxtrue+1):end,:),[model.nxtrue,np,np]);
                    vec = sym(eye(np));
                case 2
                    s2x = reshape(sx((model.nxtrue+1):end,:),[model.nxtrue,np,1]);
                    vec = model.sym.k((end-np+1):end);
            end
            for ievent = 1:nevent
                
                this.sym(ievent,:,:) = (jacobian(model.fun.sroot.sym(ievent,:),p) + jacobian(model.fun.sroot.sym(ievent,:),x(1:model.nxtrue))*sx(1:model.nxtrue,:) + jacobian(model.fun.sroot.sym(ievent,:),x(1:model.nxtrue))*sx(1:model.nxtrue,:))*vec;
                for ix = 1:model.nxtrue
                    this.sym(ievent,:,:) = this.sym(ievent,:,:) + model.fun.drootdx.sym(ievent,ix)*s2x(ix,:,:);
                end
            end
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
            % create cell array of same size
            staus = cell(1,np);
            % fill cells
            for j=1:np
                staus{j} = sprintf('stau_%i',j-1);
            end
            % transform to symbolic variable
            staus = sym(staus);
            % multiply
            this.strsym = staus;
            
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
                    dtdp = model.fun.stau.strsym; % this 1 here is correct, we explicitely do not want ievent here as the actual stau_tmp will only have dimension np
                    
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
            if (model.nxtrue < nx)
                ng_tmp = round(nx / model.nxtrue);
                this.sym = sym(zeros(np*ng_tmp,nevent));
            else
                this.sym = sym(zeros(np,nevent));
            end
            
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
            model.z2event = zeros(length(this.sym),1);
            iz = 0;
            for ievent = 1:nevent
                for jz = 1:length(model.event(ievent).z)
                    iz = iz+1;
                    model.z2event(iz) = ievent;
                end
            end
            this = makeStrSymsFull(this);
            
        case 'dzdp'
            this.sym = jacobian(model.fun.z.sym,p);
            % create cell array of same size
            this = makeStrSyms(this);
            
        case 'dzdx'
            this.sym = jacobian(model.fun.z.sym,x);
            this = makeStrSyms(this);
            
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
            this.sym = model.sym.Jy;
            % replace unify symbolic expression
            this.sym = subs(this.sym,model.sym.y,model.fun.y.strsym);
        case 'dJydy'
            this.sym = sym(zeros(model.nytrue, model.ng, model.ny));
            for iy = 1 : model.nytrue
                this.sym(iy,:,:) = jacobian(model.fun.Jy.sym(iy,:),model.fun.y.strsym);
            end
            this = makeStrSyms(this);
        case 'dJydx'
            this.sym = sym(zeros(model.nytrue, model.nxtrue, model.ng)); % Maybe nxtrue is sufficient...
            dJydy_tmp = sym(zeros(model.ng, model.ny));
            for iy = 1 : model.nytrue
                dJydy_tmp(:,:) = model.fun.dJydy.sym(iy,:,:);
                this.sym(iy,:,:) = transpose(dJydy_tmp * model.fun.dydx.strsym(:,1:model.nxtrue));
                % Transposition is necessary to have things sorted
                % correctly in gccode.m
            end
            disp('');
        case 'dJydsigma'
            this.sym = sym(zeros(model.nytrue, model.ng, length(model.sym.sigma_y)));
            for iy = 1 : model.nytrue
                this.sym(iy,:,:) = jacobian(model.fun.Jy.sym(iy,:),model.fun.sigma_y.strsym);
            end
        case 'dJydp'
            this.sym = sym(zeros(model.nytrue, model.np, model.ng));
            dJydy_tmp = sym(zeros(model.ng, model.ny));
            dJydsigma_tmp = sym(zeros(model.ng, length(model.sym.sigma_y)));
            for iy = 1 : model.nytrue
                dJydy_tmp(:,:) = model.fun.dJydy.sym(iy,:,:);
                dJydsigma_tmp(:,:) = model.fun.dJydsigma.sym(iy,:,:);
                this.sym(iy,:,:) = transpose(dJydy_tmp * model.fun.dydp.strsym ...
                    + dJydsigma_tmp * model.fun.dsigma_ydp.strsym);
                % Transposition is necessary to have things sorted
                % correctly in gccode.m
            end            
            this = makeStrSyms(this);
        case 'sJy'
            this.sym = sym(zeros(model.nytrue, model.np, model.ng));
            dJydy_tmp = sym(zeros(model.ng, model.ny));
            for iy = 1 : model.nytrue
                dJydy_tmp(:,:) = model.fun.dJydy.strsym(iy,:,:);
                this.sym(iy,:,:) = transpose(dJydy_tmp*model.fun.sy.strsym);
                % Transposition is necessary to have things sorted
                % correctly in gccode.m
            end
            this.sym = this.sym + model.fun.dJydp.strsym;
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
    if(ismember('x',this.deps))
        this.sym = mysubs(this.sym,model.sym.x,model.fun.x.sym);
    end
    this.sym = mysubs(this.sym,model.sym.p,model.fun.p.sym);
    this.sym = mysubs(this.sym,model.sym.k,model.fun.k.sym);
end

function this = makeStrSyms(this)
    this.strsym = sym(zeros(size(this.sym)));
    idx = find(this.sym);
    idx = transpose(idx(:));
    for isym = idx
        this.strsym(isym) = sym(sprintf([this.cvar '_%i'], isym-1));
    end
end

function this = makeStrSymsFull(this)
    this.strsym = sym(zeros(size(this.sym)));
    for isym = 1:numel(this.strsym)
        this.strsym(isym) = sym(sprintf([this.cvar '_%i'], isym-1));
    end
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