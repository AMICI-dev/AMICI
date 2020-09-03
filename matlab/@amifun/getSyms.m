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
                xs{j} = sprintf('var_x_%i',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(xs);
            x = this.sym;
            
        case 'dx'
            % create cell array of same size
            dxs = cell(nx,1);
            % fill cell array
            for j=1:nx
                dxs{j} = sprintf('var_dx_%i',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(dxs);
            
        case 'p'
            % create cell array of same size
            ps = cell(np,1);
            % fill cell array
            for j=1:np
                ps{j} = sprintf('var_p_%i',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(ps);
            p = this.sym;
            
        case 'k'
            % create cell array of same size
            ks = cell(nk,1);
            % fill cell array
            for j=1:nk
                ks{j} = sprintf('var_k_%i',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(ks);
            
        case 'sx'
            % create cell array of same size
            sxs = cell(nx,1);
            % fill cell array
            for j = 1:nx
                sxs{j} = sprintf('var_sx_%i', j-1);
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
                    sdx{j,i} = sprintf('var_sdx_%i', j-1);
                end
            end
            % transform into symbolic expression
            this.sym = sym(sdx);
            
        case 'xB'
            % create cell array of same size
            xBs = cell(nx,1);
            % fill cell array
            for j = 1:nx
                xBs{j} = sprintf('var_xB_%i', j-1);
            end
            % transform into symbolic expression
            this.sym = sym(xBs);
            
        case 'dxB'
            % create cell array of same size
            dxBs = cell(nx,1);
            % fill cell array
            for j = 1:nx
                dxBs{j} = sprintf('var_dxB_%i', j-1);
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
            this = makeStrSymsFull(this);
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'sigma_z'
            this.sym = model.sym.sigma_z;
            this = makeStrSymsFull(this);
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'M'
            this.sym = sym(model.sym.M);
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
                ws{iw} = sprintf('var_w_%i', iw-1);
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
            else
                model.updateRHS(tmpxdot); % update RHS anyways to generate sym_noopt for fun.xdot
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
                        this.strsym(idx_w(iw)) = sym(sprintf('var_dwdx_%i', iw-1));
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
                        this.strsym(idx_w(iw)) = sym(sprintf('var_dwdp_%i', iw-1));
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
            
            this = makeStrSymsSparse(this);
               
        
            
        case 'JDiag'
            this.sym = diag(model.fun.J.sym);
            this = makeStrSyms(this);
            
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
                dxdotdps{ix} = sprintf('tmp_dxdotdp_%i',ix-1);
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
                    this.sym=model.fun.J.strsym*sx(:,1)-model.fun.M.strsym*model.fun.sdx.sym(:,1)+model.fun.dxdotdp.strsym;
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
            
        case 'xBdot'
            if(strcmp(model.wtype,'iw'))
                syms t
                this.sym = diff(transpose(model.fun.M.sym),t)*model.fun.xB.sym + transpose(model.fun.M.sym)*model.fun.dxB.sym - transpose(model.fun.dfdx.sym)*model.fun.xB.sym;
            else
                this.sym = model.fun.JB.sym * model.fun.xB.sym;
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
            for ir = 1:nevent
                if(not(all([model.splineflag,model.minflag,model.maxflag])))
                    str = char(this.sym(ir));
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
            
        case 'drootdp'
            this.sym = jacobian(model.fun.root.sym,p);
            
        case 'drootdx'
            this.sym = jacobian(model.fun.root.sym,x);
            
        case 'drootdt'
            % noopt is important here to get derivatives right
            this.sym = diff(model.fun.root.sym,'t') + model.fun.drootdx.sym*model.fun.xdot.sym_noopt;
            
        case 'sroot'
            this.sym = model.fun.drootdp.sym + model.fun.drootdx.sym*sx;
            
        case 'srz'
            if(isfield(model.sym,'rz')) % user defined input or from augmentation
                this.sym = jacobian(model.fun.rz.sym,p) + jacobian(model.fun.rz.sym,x)*sx;
            else
                for iz = 1:length(model.z2event)
                    this.sym(iz,:) = model.fun.sroot.sym(model.z2event(iz),:);
                end
            end
          
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
            
        case 'dtaudx'
            this.sym = sym(zeros(nevent,nx));
            for ievent = 1:nevent
                this.sym(ievent,:) = - model.fun.drootdx.sym(ievent,:)/model.fun.drootdt.sym(ievent);
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
                staus{j} = sprintf('var_stau_%i',j-1);
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
            if(nx>0)
                for ievent = 1:nevent
                    this.sym(:,ievent,:) = jacobian(model.fun.deltax.sym(:,ievent),x);
                end
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
                        dxdp = sym(zeros(nx,np));
                        for ix = 1:nx
                            dxdp(ix,:) = model.fun.xdot.sym(ix)*dtdp + sx(ix,:);
                        end
                        
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
                this.sym(1:np,ievent) = -transpose(model.fun.xB.sym)*squeeze(model.fun.ddeltaxdp.sym(:,ievent,:));
                % This is just a very quick fix. Events in adjoint systems
                % have to be implemented in a way more rigorous way later
                % on... Some day...
            end
            
        case 'deltaxB'
            this.sym = sym(zeros(nx,nevent));
            for ievent = 1:nevent
                this.sym(:,ievent) = -transpose(squeeze(model.fun.ddeltaxdx.sym(:,ievent,:)))*model.fun.xB.sym;
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
            
        case 'rz'
            this.sym = sym(zeros(size(model.fun.z.sym)));
            if(isfield(model.sym,'rz'))
                this.sym = model.sym.rz;
            else
                for iz = 1:length(model.z2event)
                    this.sym(iz) = model.fun.root.sym(model.z2event(iz));
                end
            end
            this = unifySyms(this,model);
            this = makeStrSymsFull(this);
            
        case 'dzdp'
            this.sym = jacobian(model.fun.z.sym,p);
                
            for iz = 1:nz
                this.sym(iz,:) = this.sym(iz,:) + diff(model.fun.z.sym(iz),sym('t'))*model.fun.dtaudp.sym(model.z2event(iz),:);
            end
            % create cell array of same size
            this = makeStrSyms(this);
            
        case 'drzdp'
            this.sym = jacobian(model.fun.rz.sym,p);
            this = makeStrSyms(this);
            
        case 'dzdx'
            this.sym = jacobian(model.fun.z.sym,x);
            for iz = 1:nz
                this.sym(iz,:) = this.sym(iz,:)+ diff(model.fun.z.sym(iz),sym('t'))*model.fun.dtaudx.sym(model.z2event(iz),:);
            end
            this = makeStrSyms(this);
            
        case 'drzdx'
            this.sym = jacobian(model.fun.rz.sym,x);
            this = makeStrSyms(this);
            
        case 'dzdt'
            if(nz>0)
                this.sym = diff(model.fun.z.sym,'t')+jacobian(model.fun.z.sym,x(1:model.nxtrue))*model.fun.xdot.sym_noopt(1:model.nxtrue);
            else
                this.sym = sym.empty();
            end
            
        case 'sz'
            this.sym = sym(zeros(nz,np));
            tmpsym = sym(zeros(nz-model.nztrue,np));
            for iz = 1:nz
                if iz <= model.nztrue
                    this.sym(iz,:) = ...
                        + jacobian(model.fun.z.sym(iz),p) ...
                        + jacobian(model.fun.z.sym(iz),x(1:model.nxtrue))*sx(1:model.nxtrue,:) ...
                        + model.fun.dzdt.sym(iz)*model.fun.stau.sym(model.z2event(iz),:);
                else
                    % I have discovered a truly marvelous proof of these equations, which this margin is too small to
                    % contain.
                    %
                    % we need (dtau/dt)\(ddzdpdt*dtdp + dtdpddzdtdp + ddzdpdp       + dtdp*ddzdtdt*dtdp here
                    % term above:        jacx*sx+jacp   dzdt*st       jacx*sx+jacp    dzdt*st
                    % term in aug z:     drootdt        dzdp          dzdp            drootdt
                    % augmented is always drootdt\dzdp or linear combinations
                    % we cannot correctly compute ddzdpdt, but only ddzdtdp so we omit the drootdt part of jacx*sx-jacp
                    % symmetrise it and add it later on.
                    % also the dzdp part contains second order sensitivities and the drootdt part does not (this leads
                    % to 1:model.nxtrue indexing)
                    
                    
                    if(model.o2flag==1)
                        tmpsym(iz-model.nztrue,:) = jacobian(1/model.fun.drootdt.sym(model.z2event(iz)),p)*model.fun.z.sym(iz)*model.fun.drootdt.sym(model.z2event(iz)) ...
                            + jacobian(1/model.fun.drootdt.sym(model.z2event(iz)),x(1:model.nxtrue))*sx(1:model.nxtrue,:)*model.fun.z.sym(iz)*model.fun.drootdt.sym(model.z2event(iz));
                    else
                       error('sz for directional second order sensis was never implemented and I do not know how to, you are on your own here.'); 
                    end
                    
                    this.sym(iz,:) = ...
                        + jacobian(model.fun.z.sym(iz)*model.fun.drootdt.sym(model.z2event(iz)),p)/model.fun.drootdt.sym(model.z2event(iz)) ...
                        + jacobian(model.fun.z.sym(iz)*model.fun.drootdt.sym(model.z2event(iz)),x)*sx/model.fun.drootdt.sym(model.z2event(iz)) ...
                        + model.fun.dzdt.sym(iz)*model.fun.stau.sym(model.z2event(iz),:);
                end
            end
            if(model.nz>0)
                if(model.o2flag==1)
                    % THIS IS WHAT YOU NEED TO FIX FOR model.o2flag == 2, transpose does not do the deal :(
                    % we need to pay attention here, some sensitivities are encoded as sx and some as augmented x.
                    % the sx ones are not transposable (as the parameter is encoded in the matrix column) so we need to
                    % convert all of them to augmented x.
                    for ip = 1:np
                        tmpsym(:,ip) = subs(tmpsym(:,ip),sx(1:model.nxtrue,1),x((model.nxtrue*ip+1):(model.nxtrue*(ip+1))));
                    end
                    this.sym(model.nztrue+1:end,:) = this.sym(model.nztrue+1:end,:) + tmpsym + transpose(tmpsym);
                    % you might not believe this, but this matrix should (and hopefully will) actually be symmetric ;)
                end
            end
            
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
            
        case 'JBand'
            %do nothing
        case 'JBandB'
            %do nothing
        case 'JSparse'
            %do nothing
            
        case 'Jy'
            this.sym = model.sym.Jy;
            % replace unify symbolic expression
            this = unifySyms(this,model);
            tmp = arrayfun(@(iy)sym(sprintf('y_%i',iy-1)),1:ny,'UniformOutput',false);
            this.sym = mysubs(this.sym,[tmp{:}],model.fun.y.strsym);
            tmp = arrayfun(@(iy)sym(sprintf('sigma_y_%i',iy-1)),1:ny,'UniformOutput',false);
            this.sym = mysubs(this.sym,[tmp{:}],model.fun.sigma_y.strsym);
        case 'dJydy'
            this.sym = sym(zeros(model.nytrue, model.ng, model.ny));
            for iyt = 1 : model.nytrue
                this.sym(iyt,:,:) = jacobian(model.fun.Jy.sym(iyt,:),model.fun.y.strsym);
            end
        case 'dJydsigma'
            this.sym = sym(zeros(model.nytrue, model.ng, model.ny));
            for iyt = 1 : model.nytrue
                this.sym(iyt,:,:) = jacobian(model.fun.Jy.sym(iyt,:),model.fun.sigma_y.strsym);
            end
        case 'Jz'
            this.sym = model.sym.Jz;
            this = unifySyms(this,model);
            z = arrayfun(@(iz)sym(sprintf('z_%i',iz-1)),1:nz,'UniformOutput',false);
            this.sym = mysubs(this.sym,[z{:}],model.fun.z.strsym);
            tmp = arrayfun(@(iz)sym(sprintf('sigma_z_%i',iz-1)),1:nz,'UniformOutput',false);
            this.sym = mysubs(this.sym,[tmp{:}],model.fun.sigma_z.strsym);
        case 'Jrz'
            this.sym = model.sym.Jrz;
            this = unifySyms(this,model);
            tmp = arrayfun(@(iz)sym(sprintf('rz_%i',iz-1)),1:nz,'UniformOutput',false);
            this.sym = mysubs(this.sym,[tmp{:}],model.fun.z.strsym);
            tmp = arrayfun(@(iz)sym(sprintf('sigma_z_%i',iz-1)),1:nz,'UniformOutput',false);
            this.sym = mysubs(this.sym,[tmp{:}],model.fun.sigma_z.strsym);
        case 'dJzdz'
            this.sym = sym(zeros(model.nztrue, model.ng, model.nz));
            for iz = 1 : model.nztrue
                this.sym(iz,:,:) = jacobian(model.fun.Jz.sym(iz,:),model.fun.z.strsym);
            end
            this = makeStrSyms(this);
        case 'dJrzdz'
            this.sym = sym(zeros(model.nztrue, model.ng, model.nz));
            for iz = 1 : model.nztrue
                this.sym(iz,:,:) = jacobian(model.fun.Jrz.sym(iz,:),model.fun.rz.strsym);
            end
            this = makeStrSyms(this);
        case 'dJzdsigma'
            this.sym = sym(zeros(model.nztrue, model.ng, model.nz));
            for iz = 1 : model.nztrue
                this.sym(iz,:,:) = jacobian(model.fun.Jz.sym(iz,:),model.fun.sigma_z.strsym);
            end
        case 'dJrzdsigma'
            this.sym = sym(zeros(model.nztrue, model.ng, model.nz));
            for iz = 1 : model.nztrue
                this.sym(iz,:,:) = jacobian(model.fun.Jrz.sym(iz,:),model.fun.sigma_z.strsym);
            end
            
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

function this = makeStrSymsSparse(this)
    this.strsym = sym(zeros(size(this.sym)));
    idx = find(this.sym);
    for isym = 1:length(idx)
        this.strsym(idx(isym)) = sym(sprintf([this.cvar '_%i'], isym-1));
    end
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