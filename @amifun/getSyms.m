function this = getSyms(this,model)
    % getSyms computes the symbolic expression for the requested function
    %
    % Return values:
    %  this: updated function definition object @type amifun
    
    % store often used variables for ease of notation, dependencies should
    % ensure that these variables are defined
    
    persistent x p sx
    
    nx = model.nx;
    nevent = model.ne;
    np = model.np;
    nk = model.nk;
    
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
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'sigma_z'
            this.sym = model.sym.sigma_z;
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
        case 'M'
            this.sym = model.sym.M;
            % replace unify symbolic expression
            this = unifySyms(this,model);
            
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
            
        case 'dfdx'
            this.sym=jacobian(model.fun.xdot.sym,x);
            
        case 'J'
            if(strcmp(model.wtype,'iw'))
                syms cj
                this.sym = model.fun.dfdx.sym - cj*model.fun.M.sym;
            else
                this.sym = jacobian(model.fun.xdot.sym,x);
            end
            
            %% 
            % build short strings for reuse of jacobian
            
            % find nonzero entries
            idx = find(double(this.sym~=0));
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
            this.sym=jacobian(model.fun.xdot.sym,p);
            
            %%
            % build short strings for reuse of dxdotdp
            
            % find nonzero rows, only completely nonzero rows matter as we
            % will have parameter dependent implementation later on
            % anyways.
            idx = find(any(double(this.sym~=0),2));
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
            
        case 'dydp'
            this.sym=jacobian(model.fun.y.sym,p);
            
        case 'sy'
            this.sym=model.fun.dydp.sym + model.fun.dydx.sym*model.fun.sx.sym ;
            
        case 'Jv'
            % create cell array of same size
            vs = cell(nx,1);
            % fill cells
            for j=1:nx
                vs{j} = sprintf('v[%i]',j-1);
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
                vs{j} = sprintf('v[%i]',j-1);
            end
            % transform to symbolic variable
            vs = sym(vs);
            % multiply
            this.sym = -transpose(model.fun.J.sym)*vs;
            
        case 'xBdot'
            if(strcmp(model.wtype,'iw'))
                syms t
                this.sym = diff(transpose(model.fun.M.syms),t)*model.fun.xB.syms + transpose(model.fun.M.syms)*model.fun.dxB.syms - transpose(model.fun.dfdx.syms)*model.fun.xB.syms;
            else
                this.sym = -transpose(model.fun.J.sym)*model.fun.xB.sym;
            end
            
        case 'qBdot'
            this.sym = -transpose(model.fun.xB.sym)*model.fun.dxdotdp.strsym;
            
        case 'dsigma_ydp'
            this.sym = jacobian(model.fun.sigma_y.sym,p);
            
        case 'dsigma_zdp'
            this.sym = jacobian(model.fun.sigma_z.sym,p);
            
        case 'root'
            this.sym = [model.event.trigger];
            this = unifySyms(this,model);
            
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
            this.sym = [model.event.bolus];
            this = unifySyms(this,model);
            
        case 'deltaxdot'
            this.sym = [model.event.deltaxdot];
            this = unifySyms(this,model);
            
        case 'ddeltaxdp'
            this.sym = jacobian(model.fun.deltax.sym,p);
            
        case 'ddeltaxdx'
            this.sym = jacobian(model.fun.deltax.sym,x);
            
        case 'ddeltaxdt'
            this.sym = diff(model.fun.deltax.sym,'t');
            
        case 'deltasx'
            
            if(nevent>0)
                for ievent = 1:nevent
                    
                    % dtdp  = (1/drdt)*drdp
                    dtdp = model.fun.stau.sym(ievent,:);
                    
                    % dxdp  = dx/dt*dt/dp + dx/dp
                    dxdp = model.fun.xdot.sym*dtdp + sx;
                    
                    this.sym(:,:,ievent)= ...
                        + model.fun.ddeltaxdx.sym*dxdp ...
                        + model.fun.ddeltaxdt.sym*dtdp ...
                        + model.fun.ddeltaxdp.sym;
                end
            end
            
        case 'deltaqB'
            this.sym = transpose(model.fun.xB.sym)*model.fun.ddeltaxdp.sym;
                    
        case 'deltaxB'
            this.sym = model.fun.ddeltaxdx.sym*model.fun.xB.sym;
            
        case 'z'
            this.sym = [model.event.z];
            this = unifySyms(this,model);
            
        case 'dzdp'
            this.sym = jacobian(model.fun.z.sym,p);
            
        case 'dzdx'
            this.sym = jacobian(model.fun.z.sym,x);
            
        case 'dzdt'
            this.sym = diff(model.fun.z.sym,'t');
            
        case 'sz'
            this.sym = model.fun.dzdp.sym + model.fun.dzdx.sym*sx + model.fun.dzdt.sym*model.fun.stau.sym;
            
        case 'JBand'
            %do nothing
        case 'JBandB'
            %do nothing    
        case 'JSparse'
            %do nothing 
        case 'JSparseB'
            %do nothing 
            
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