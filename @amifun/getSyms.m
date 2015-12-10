function this = getSyms(this,model)
    % getfunargs returns a string which contains all input arguments for the function fun.
    
    fprintf([funstr ' | '])
    switch(this.funstr)
        case 'x'
            % create cell array of same size
            xs = cell(nx,1);
            % fill cell array
            for j=1:nx
                xs{j} = sprintf('x[%i]',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(xs);
            
        case 'dx'
            % create cell array of same size
            dxs = cell(nx,1);
            % fill cell array
            for j=1:nx
                dxs{j} = sprintf('dx[%i]',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(dxs);
            
        case 'p'
            % create cell array of same size
            ps = cell(np,1);
            % fill cell array
            for j=1:np
                ps{j} = sprintf('p[%i]',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(ps);
            
        case 'k'
            % create cell array of same size
            ks = cell(nk,1);
            % fill cell array
            for j=1:nk
                ks{j} = sprintf('k[%i]',j-1);
            end
            % transform into symbolic expression
            this.sym = sym(ks);
            
        case 'sx'
            % create cell array of same size
            sx = cell(nx,np);
            % fill cell array
            for j = 1:nx
                for i = 1:np
                    sx{j,i} = sprintf('sx[%i]', j-1);
                end
            end
            % transform into symbolic expression
            this.sym = sym(sx);
            
        case 'sdx'
            % create cell array of same size
            sdx = cell(nx,np);
            % fill cell array
            for j = 1:nx
                for i = 1:np
                    sdx{j,i} = sprintf('sdx[%i]', j-1);
                end
            end
            % transform into symbolic expression
            this.sym = sym(sdx);
            
        case 'xB'
            % create cell array of same size
            xBs = cell(nx,1);
            % fill cell array
            for j = 1:nx
                xBs{j} = sprintf('xB[%i]', j-1);
            end
            % transform into symbolic expression
            this.sym = sym(xBs);
            
        case 'dxB'
            % create cell array of same size
            dxBs = cell(nx,1);
            % fill cell array
            for j = 1:nx
                dxBs{j} = sprintf('dxB[%i]', j-1);
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
            
        case 'sigma_t'
            this.sym = model.sym.sigma_t;
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
            
            if(strcmp(this.wtype,'iw'))
                if(size(this.sym,2)>size(this.sym,1))
                    this.sym = -transpose(model.funs.M.sym*model.funs.dx.sym)+this.sym;
                else
                    this.sym = -model.funs.M.sym*model.funs.dx.sym+this.sym;
                end
                
                this.id = sum(model.funs.M.sym,2)~=0;
            else
                this.id = zeros(this.nx,1);
            end
            
        case 'dfdx'
            this.sym=jacobian(model.funs.xdot.sym,model.funs.x.sym);
            
        case 'J'
            if(strcmp(this.wtype,'iw'))
                syms cj
                this.sym = model.funs.dfdx.sym - cj*model.funs.M.sym;
            else
                this.sym.J = jacobian(this.sym.xdot,this.strsym.xs);
            end
            
            %% 
            % build short strings for reuse of jacobian
            
            % find nonzero entries
            idx = find(double(this.sym~=0));
            % create cell array of same size
            Js = cell(length(idx),1);
            % fill cells with strings
            for iJ = 1:length(idx)
                Js{iJ} = sprintf('tmp_J[%i]',iJ-1);
            end
            % create full symbolic matrix
            this.strsym = sym(zeros(nx,nx));
            % fill nonzero entries
            this.strsym(idx) = sym(Js);
            
        case 'JB'
            this.sym.JB=-transpose(this.sym.J);
            
        case 'dxdotdp'
            this.sym=jacobian(model.funs.xdot.sym,model.funs.p.sym);
            
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
                dxdotdps{ix} = sprintf('tmp_dxdotdp[%i+ip*nx]',ix-1);
            end
            % create full symbolic array
            this.strsym = sym(zeros(nx,1));
            % fill nonzero entries
            this.strsym(idx) = sym(dxdotdps);
            
        case 'sx0'
            this.sym=jacobian(model.funs.x0.sym,model.funs.p.sym);
            
        case 'sdx0'
            this.sym.sdx0=jacobian(model.funs.dx0.sym,model.funs.p.sym);
            
        case 'sxdot'
            if(strcmp(this.wtype,'iw'))
                this.sym=model.funs.dfdx.sym*model.funs.sx.sym-model.funs.M.sym*model.funs.sdx.sym+model.funs.dxdotdp.sym;
            else
                this.sym=model.funs.J.strsym*model.funs.sx.sym(:,1)+model.funs.dxdotdp.sym;
            end
            
        case 'dydx'
            this.sym=jacobian(model.funs.y.sym,model.funs.x.sym);
            
        case 'dydp'
            this.sym=jacobian(model.funs.y.sym,model.funs.p.sym);
            
        case 'sy'
            this.sym=model.funs.dydp.sym + this.sym.dydx*this.strsym.sx ;
            
            % events
        case 'drootdx'
            if(this.nxtrue > 0)
                this.sym.drootdx = jacobian(this.sym.rfun,this.strsym.xs(1:this.nxtrue));
            else
                this.sym.drootdx=jacobian(this.sym.rfun,this.strsym.xs);
            end
            
        case 'drootdt'
            if(this.nxtrue > 0)
                this.sym.drootdt = diff(this.sym.rfun,sym('t')) ...
                    + this.sym.drootdx*this.sym.xdot(1:this.nxtrue);
            else
                this.sym.drootdt = diff(this.sym.rfun,sym('t')) + this.sym.drootdx*this.sym.xdot;
            end
            
        case 'drootpdp'
            this.sym.drootpdp = jacobian(this.sym.rfun,this.strsym.ps);
            
        case 'drootdp'
            if(this.nxtrue > 0)
                this.sym.drootdp = this.sym.drootpdp + this.sym.drootdx*this.sym.dxdp;
            else
                this.sym.drootdp = this.sym.drootpdp + this.sym.drootdx*this.strsym.sx;
            end
            
        case 'sroot'
            this.sym.sroot = sym.zeros(nr,np);
            if(this.nr>0)
                for ir = 1:this.nr
                    this.sym.sroot(ir,:) = -(this.sym.drootdt(ir,:))\(this.sym.drootdp(ir,:));
                end
            end
            
        case 's2root'
            this.sym.s2root =  sym(zeros(nr,np,np));
            % ddtaudpdp = -dgdt\(dtaudp*ddgdtdt*dtaudp + 2*ddgdpdt*dtaudp
            % + ddgdpdp)
            if(this.nr>0 && this.nxtrue > 0)
                for ir = 1:nr
                    this.sym.s2root(ir,:,:) = -1/(this.sym.drootdt(ir))*(...
                        transpose(this.sym.sroot(ir,:))*this.sym.ddrootdtdt(ir)*this.sym.sroot(ir,:) ...
                        + squeeze(this.sym.ddrootdtdp(ir,:,:))*this.sym.sroot(ir,:) ...
                        + transpose(squeeze(this.sym.ddrootdtdp(ir,:,:))*this.sym.sroot(ir,:)) ...
                        + squeeze(this.sym.s2rootval(ir,:,:)));
                end
            end
            
        case 'srootval'
            if(this.nr>0)
                this.sym.srootval = this.sym.drootdp;
            else
                this.sym.srootval = sym(zeros(0,np));
            end
            
        case 's2rootval'
            if(this.nr>0 && this.nxtrue > 0)
                this.sym.s2rootval = this.sym.ddrootdpdp;
            else
                this.sym.s2rootval =  sym(zeros(nr,np,np));
            end
            
        case 'dxdp'
            if(this.nr>0 && this.nxtrue > 0)
                this.sym.dxdp = reshape(this.strsym.xs((this.nxtrue+1):end),this.nxtrue,np);
            else
                this.sym.dxdp = sym(zeros(this.nxtrue,np));
            end
        case 'ddrootdxdx'
            this.sym.ddrootdxdx = sym(zeros(nr,this.nxtrue,this.nxtrue));
            if(this.nxtrue>0)
                for ir = 1:nr
                    this.sym.ddrootdxdx(ir,:,:) = hessian(this.sym.rfun(ir),this.strsym.xs(1:this.nxtrue));
                end
            end
            
        case 'ddrootdxdt'
            this.sym.ddrootdxdt = sym(zeros(nr,this.nxtrue));
            % ddgdxdt = ddgdxdt + ddgdxdx*dxdt
            if(this.nxtrue>0)
                for ir = 1:nr
                    this.sym.ddrootdxdt(ir,:) = diff(this.sym.drootdx(ir,:),sym('t')) ...
                        + jacobian(this.sym.drootdx*this.sym.xdot(1:this.nxtrue),this.strsym.xs(1:this.nxtrue)) ...
                        + transpose(squeeze(this.sym.ddrootdxdx(ir,:,:))*this.sym.xdot(1:this.nxtrue));
                end
            end
            
        case 'ddrootdtdt'
            % ddgdtdt = pddgdtpdt + ddgdxdt*dxdt
            if(this.nxtrue>0)
                this.sym.ddrootdtdt = diff(this.sym.drootdt,sym('t')) ...
                    + this.sym.ddrootdxdt*this.sym.xdot(1:this.nxtrue);
            else
                this.sym.ddrootdtdt = sym(zeros(size(this.sym.rfun)));
            end
            
        case 'ddxdtdp'
            % ddxdtdp = pddxdtpdp + ddxdtdx*dxdp
            if(this.nxtrue>0)
                this.sym.ddxdtdp = this.sym.dxdotdp(1:this.nxtrue,:) ...
                    + this.sym.dxdotdx*this.sym.dxdp;
            else
                this.sym.ddxdtdp = sym(zeros(nx,np));
            end
            
        case 'dxdotdx'
            this.sym.dxdotdx = jacobian(this.sym.xdot(1:this.nxtrue),this.strsym.xs(1:this.nxtrue));
            
        case 'ddrootdpdx'
            this.sym.ddrootdpdx = sym(zeros(nr,np,this.nxtrue));
            if(this.nxtrue>0)
                for ir = 1:nr
                    this.sym.ddrootdpdx(ir,:,:) = jacobian(this.sym.drootdp(ir,:),this.strsym.xs(this.rt(1:this.nxtrue)));
                end
            end
            
        case 'ddrootdpdt'
            % ddgdpdt =  pddgdppdt + ddgdpdx*dxdt
            this.sym.ddrootdpdt = sym(zeros(nr,np,1));
            if(this.nxtrue>0)
                for ir = 1:nr
                    this.sym.ddrootdpdt(ir,:,:) = diff(this.sym.drootdp,sym('t')) ...
                        + transpose(permute(this.sym.ddrootdpdx(ir,:,:),[2,3,1])*this.sym.xdot(1:this.nxtrue));
                end
            end
            
        case 'ddrootdtdp'
            % ddgdpdt =  pddgdtpdp + ddgdtdx*dxdp
            this.sym.ddrootdtdp = sym(zeros(nr,1,np));
            if(this.nxtrue>0)
                for ir = 1:nr
                    this.sym.ddrootdtdp(ir,:,:) = jacobian(this.sym.drootdt,this.strsym.ps) ...
                        + this.sym.ddrootdxdt*this.sym.dxdp;
                end
            end
            
        case 'drootdx_ddxdpdp'
            % dgdx*ddxdpdp
            this.sym.drootdx_ddxdpdp = sym(zeros(nr,np,np));
            if(this.nxtrue>0)
                ddxdpdp = reshape(this.strsym.sx((this.nxtrue+1):end,:),this.nxtrue,np,np);
                for ir = 1:nr
                    for ix = 1:this.nxtrue
                        this.sym.drootdx_ddxdpdp(ir,:,:) = this.sym.drootdx_ddxdpdp(ir,:,:) + this.sym.drootdx(ir,ix)*ddxdpdp(ix,:,:);
                    end
                end
            end
            
        case 'ddrootdpdp'
            % ddgdpdp = pddgdppdp + ddgdpdx*dxdp + dgdx*dxdp
            this.sym.ddrootdpdp = sym(zeros(nr,np,np));
            if(this.nxtrue>0)
                for ir = 1:nr
                    this.sym.ddrootdpdp(ir,:,:) = jacobian(this.sym.drootdp,this.strsym.ps) ...
                        + permute(this.sym.ddrootdpdx(ir,:,:),[2,3,1])*this.sym.dxdp ...
                        + squeeze(this.sym.drootdx_ddxdpdp(ir,:,:));
                end
            end
            
            
        case 'dtdp'
            this.sym.dtdp = sym(zeros(nr,np));
            if(nr>0)
                for ir = 1:nr
                    this.sym.dtdp(ir,:) = -(this.sym.drootdt(ir,:))\(this.sym.drootpdp(ir,:));
                end
            end
            
        case 'dtdx'
            this.sym.dtdx = sym(zeros(nr,nx));
            if(nr>0)
                for ir = 1:nr
                    this.sym.dtdx(ir,:) = -(this.sym.drootdt(ir,:))\(this.sym.drootdx(ir,:));
                end
            end
            
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
            this.sym = model.funs.J.sym*vs;
            
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
            this.sym = -transpose(model.funs.J.sym)*vs;
            
        case 'xBdot'
            if(strcmp(this.wtype,'iw'))
                syms t
                this.sym = diff(transpose(model.funs.M.syms),t)*model.funs.xB.syms + transpose(model.funs.M.syms)*model.funs.dxB.syms - transpose(model.funs.dfdx.syms)*model.funs.xB.syms;
            else
                this.sym = -transpose(model.funs.J.syms)*model.funs.xB.syms;
            end
            
        case 'qBdot'
            this.sym = -transpose(model.funs.xB.syms)*model.funs.dxdotdp.strsyms;
            
        case 'dsigma_ydp'
            this.sym = jacobian(model.funs.sigma_y.syms,model.funs.sigma_p.syms);
            
        case 'dsigma_tdp'
            this.sym = jacobian(model.funs.sigma_t.syms,model.funs.sigma_p.syms);
            
        case 'rdisc'
            % bogus variable
            syms polydirac
            % discontinuities
            ndisc = 0;
            for ix = 1:length(this.sym.xdot)
                tmp_str = char(this.sym.xdot(ix));
                idx_start = [strfind(tmp_str,'heaviside') + 9,strfind(tmp_str,'dirac') + 5]; % find
                if(~isempty(idx_start))
                    for iocc = 1:length(idx_start) % loop over occurances
                        brl = 1; % init bracket level
                        str_idx = idx_start(iocc); % init string index
                        % this code should be improved at some point
                        while brl >= 1 % break as soon as initial bracket is closed
                            str_idx = str_idx + 1;
                            if(strcmp(tmp_str(str_idx),'(')) % increase bracket level
                                brl = brl + 1;
                            end
                            if(strcmp(tmp_str(str_idx),')')) % decrease bracket level
                                brl = brl - 1;
                            end
                        end
                        idx_end = str_idx;
                        arg = tmp_str((idx_start(iocc)+1):(idx_end-1));
                        if(ndisc>0)
                            ediscf = zeros(size(this.sym.rdisc));
                            for id = 1:length(this.sym.rdisc)
                                ediscf(id) = isequaln(abs(this.sym.rdisc(id)),abs(sym(arg)));
                            end
                            if(~any(ediscf))
                                if(isempty(strfind(arg,','))) % no higher order derivatives
                                    ndisc = ndisc + 1;
                                    this.sym.rdisc(ndisc) = sym(arg);
                                end
                            else
                            end
                        else
                            if(isempty(strfind(arg,','))) % no higher order derivatives
                                ndisc = ndisc + 1;
                                this.sym.rdisc(ndisc) = sym(arg);
                            end
                        end
                    end
                end
            end
            this.sym.f = sym(zeros(nx,ndisc));
            for ix = 1:nx
                symchar = char(this.sym.xdot(ix));
                if(strfind(symchar,'dirac'))
                    for idisc = 1:length(this.sym.rdisc)
                        discchar = char(this.sym.rdisc(idisc));
                        str_arg_d = ['dirac(' discchar ')' ];
                        str_arg_d1 = ['dirac(1, ' discchar ')' ];
                        str_arg_1d = ['dirac(' discchar ', 1)' ];
                        str_arg_d2 = ['dirac(2, ' discchar ')' ];
                        str_arg_2d = ['dirac(' discchar ', 2)' ];
                        % f
                        symchar = strrep(symchar,str_arg_d,'polydirac');
                        symchar = strrep(symchar,str_arg_d1,'polydirac');
                        symchar = strrep(symchar,str_arg_1d,'polydirac');
                        symchar = strrep(symchar,str_arg_d2,'polydirac');
                        symchar = strrep(symchar,str_arg_2d,'polydirac');
                    end
                    cfp = coeffs(sym(symchar),polydirac);
                    this.sym.f(ix) = cfp(1);
                else
                    this.sym.f(ix) = this.sym.xdot(ix);
                end
            end
            for idisc = 1:length(this.sym.rdisc)
                this.sym.gdisc(ndisc) = diff(this.sym.rdisc(idisc),'t') + jacobian(this.sym.rdisc(idisc),this.strsym.xs)*this.sym.f;
            end
            this.ndisc = ndisc;
            if(ndisc == 0)
                this.sym.rdisc = [];
            end
            
        case 'deltadisc'
            syms polydirac
            
            this.sym.deltadisc = sym(zeros(nx,ndisc));
            this.sym.v = sym(zeros(nx,ndisc));
            this.sym.w = sym(zeros(nx,ndisc));
            this.sym.z = sym(zeros(nx,ndisc));
            
            % convert deltas in xdot
            if(ndisc>0)
                for ix = 1:nx
                    symchar = char(this.sym.xdot(ix));
                    if(strfind(symchar,'dirac'))
                        for idisc = 1:ndisc
                            discchar = char(this.sym.rdisc(idisc));
                            str_arg_d = ['dirac(' discchar ')' ];
                            str_arg_d1 = ['dirac(1, ' discchar ')' ];
                            str_arg_1d = ['dirac(' discchar ', 1)' ];
                            str_arg_d2 = ['dirac(2, ' discchar ')' ];
                            str_arg_2d = ['dirac(' discchar ', 2)' ];
                            % find the corresponding discontinuity
                            % xdot_i = f_i(x,t) + v_i*delta(gdisc) +
                            % w_i*delta(1)(gdisc) + z_i*delta(2)(gdisc);
                            % v
                            tmpstr = char(this.sym.xdot(ix));
                            tmpstr = strrep(tmpstr,str_arg_d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.v(ix,idisc) = cfp(2)*1/abs(this.sym.gdisc(idisc));
                                this.sym.vcoeff(ix,idisc) = cfp(2);
                            else
                                this.sym.v(ix,idisc) = 0;
                            end
                            % w
                            tmpstr = char(this.sym.xdot(ix));
                            tmpstr = strrep(tmpstr,str_arg_d1,'polydirac');
                            tmpstr = strrep(tmpstr,str_arg_1d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.w(ix,idisc) = cfp(2)*this.sym.gdisc(idisc)/abs(this.sym.gdisc(idisc));
                            else
                                this.sym.w(ix,idisc) = 0;
                            end
                            % z
                            tmpstr = char(this.sym.xdot(ix));
                            tmpstr = strrep(tmpstr,str_arg_d2,'polydirac');
                            tmpstr = strrep(tmpstr,str_arg_2d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.z(ix,idisc) = cfp(2)*this.sym.gdisc(idisc)^2/abs(this.sym.gdisc(idisc));
                            else
                                this.sym.z(ix,idisc) = 0;
                            end
                        end
                    end
                end
                
                % compute deltadisc
                for idisc = 1:ndisc
                    nonzero_z = any(this.sym.z(:,idisc)~=0);
                    if(nonzero_z)
                        M = jacobian(this.sym.f(:,idisc),this.strsym.xs)*this.sym.z(:,idisc) ...
                            + this.sym.w(:,idisc) ...
                            - diff(this.sym.z(:,idisc),'t') - jacobian(this.sym.z(:,idisc),this.strsym.xs)*this.sym.f(:,idisc);
                    else
                        M = this.sym.w(:,idisc);
                    end
                    nonzero_M = any(M~=0);
                    for ix = 1:nx
                        % delta_x_i =
                        % sum_j(df_idx_j*(sum_k(df_idx_k*z_k) +
                        % w_j - dot{z}_j) + v_i - dot{w}_i +
                        % ddot{z}_i
                        
                        if(this.sym.v(ix,idisc)~=0)
                            this.sym.deltadisc(ix,idisc) = this.sym.deltadisc(ix,idisc) ...
                                +this.sym.v(ix,idisc);
                        end
                        
                        if(this.sym.w(ix,idisc)~=0)
                            this.sym.deltadisc(ix,idisc) = this.sym.deltadisc(ix,idisc) ...
                                - diff(this.sym.w(ix,idisc),'t') - jacobian(this.sym.w(ix,idisc),this.strsym.xs)*this.sym.f(:,idisc);
                        end
                        
                        if(this.sym.z(ix,idisc)~=0)
                            this.sym.deltadisc(ix,idisc) = this.sym.deltadisc(ix,idisc) ...
                                + diff(this.sym.z(ix,idisc),'t',2) + diff(jacobian(this.sym.z(ix,idisc),this.strsym.xs)*this.sym.f(:,idisc),'t') + ...
                                jacobian(jacobian(this.sym.z(ix,idisc),this.strsym.xs)*this.sym.xdot,this.strsym.xs)*this.sym.f(:,idisc);
                        end
                        
                        if(nonzero_M)
                            this.sym.deltadisc(ix,idisc) = this.sym.deltadisc(ix,idisc) ...
                                +jacobian(this.sym.f(ix,idisc),this.strsym.xs)*M;
                        end
                        
                    end
                end
            end
            
        case 'ideltadisc'
            syms polydirac
            this.sym.ideltadisc = sym(zeros(np,ndisc));
            this.sym.if = sym(zeros(np,ndisc));
            this.sym.iv = sym(zeros(np,ndisc));
            this.sym.iw = sym(zeros(np,ndisc));
            this.sym.iz = sym(zeros(np,ndisc));
            
            % convert deltas in qBdot
            if(ndisc>0)
                for ip = 1:np
                    symchar = char(this.sym.eqBdot(ip));
                    if(strfind(symchar,'dirac'))
                        for idisc = 1:ndisc
                            discchar = char(this.sym.rdisc(idisc));
                            str_arg_d = ['dirac(' discchar ')' ];
                            str_arg_d1 = ['dirac(1, ' discchar ')' ];
                            str_arg_1d = ['dirac(' discchar ', 1)' ];
                            str_arg_d2 = ['dirac(2, ' discchar ')' ];
                            str_arg_2d = ['dirac(' discchar ', 2)' ];
                            % find the corresponding discontinuity
                            % xdot_i = f_i(x,t) + h_i*heaviside(gdisc) + v_i*delta(gdisc) +
                            % w_i*delta(1)(gdisc) + z_i*delta(2)(gdisc);
                            % v
                            tmpstr = char(this.sym.eqBdot(ip));
                            tmpstr = strrep(tmpstr,str_arg_d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.iv(ip,idisc) = cfp(2)*1/abs(this.sym.gdisc(idisc));
                            else
                                this.sym.iv(ip,idisc) = 0;
                            end
                            % w
                            tmpstr = char(this.sym.eqBdot(ip));
                            tmpstr = strrep(tmpstr,str_arg_d1,'polydirac');
                            tmpstr = strrep(tmpstr,str_arg_1d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.iw(ip,idisc) = cfp(2)*this.sym.gdisc(idisc)/abs(this.sym.gdisc(idisc));
                            else
                                this.sym.iw(ip,idisc) = 0;
                            end
                            % z
                            tmpstr = char(this.sym.eqBdot(ip));
                            tmpstr = strrep(tmpstr,str_arg_d2,'polydirac');
                            tmpstr = strrep(tmpstr,str_arg_2d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.iz(ip,idisc) = cfp(2)*this.sym.gdisc(idisc)^2/abs(this.sym.gdisc(idisc));
                            else
                                this.sym.iz(ip,idisc) = 0;
                            end
                            % f
                            symchar = strrep(symchar,str_arg_d,'polydirac');
                            symchar = strrep(symchar,str_arg_d1,'polydirac');
                            symchar = strrep(symchar,str_arg_1d,'polydirac');
                            symchar = strrep(symchar,str_arg_d2,'polydirac');
                            symchar = strrep(symchar,str_arg_2d,'polydirac');
                        end
                        cfp = coeffs(sym(symchar),polydirac);
                        this.sym.if(ip) = cfp(1);
                    else
                        this.sym.if(ip,:) = repmat(this.sym.eqBdot(ip),[1,ndisc]);
                    end
                end
                % compute ideltadisc
                for idisc = 1:ndisc
                    M = (jacobian(this.sym.f(:,idisc),this.strsym.xs)*this.sym.z(:,idisc) ...
                        + this.sym.w(:,idisc) ...
                        - diff(this.sym.z(:,idisc),'t') ...
                        - jacobian(this.sym.z(:,idisc),this.strsym.xs)*this.sym.f(:,idisc));
                    nonzero_M = any(M~=0);
                    for ip = 1:np
                        % delta_x_i =
                        % sum_j(df_idx_j*(sum_k(df_idx_k*z_k) +
                        % w_j - dot{z}_j) + v_i - dot{w}_i +
                        % ddot{z}_i
                        if(this.sym.iz(ip,idisc)~=0)
                            dizdx = jacobian(this.sym.iz(ip,idisc),this.strsym.xs);
                            dizdsxdot = jacobian(this.sym.iz(ip,idisc),this.strsym.xBs);
                            
                            this.sym.ideltadisc(ip,idisc) = this.sym.ideltadisc(ip,idisc) ...
                                + diff(this.sym.iz(ip,idisc),'t',2) ...
                                + diff(dizdx*this.sym.f(:,idisc),'t') ...
                                + diff(dizdsxdot*this.sym.xBdot,'t') ...
                                + jacobian(dizdx*this.sym.f(:,idisc),this.strsym.xs)*this.sym.f(:,idisc) ...
                                + jacobian(dizdx*this.sym.f(:,idisc),this.strsym.xBs)*this.sym.xBdot ...
                                + jacobian(dizdsxdot*this.sym.xBdot,this.strsym.xs)*this.sym.f(:,idisc) ...
                                + jacobian(dizdsxdot*this.sym.xBdot,this.strsym.xBs)*this.sym.xBdot;
                        end
                        
                        if(this.sym.iv(ip,idisc)~=0)
                            this.sym.ideltadisc(ip,idisc) = this.sym.ideltadisc(ip,idisc) ...
                                + this.sym.iv(ip,idisc);
                        end
                        
                        if(this.sym.iw(ip,idisc)~=0)
                            this.sym.ideltadisc(ip,idisc) = this.sym.ideltadisc(ip,idisc) ...
                                - diff(this.sym.iw(ip,idisc),'t') - jacobian(this.sym.iw(ip,idisc),this.strsym.xs)*this.sym.f(:,idisc) ...
                                - jacobian(this.sym.iw(ip,idisc),this.strsym.xBs)*this.sym.xBdot;
                            
                        end
                        
                        if(nonzero_M)
                            this.sym.ideltadisc(ip,idisc) = this.sym.ideltadisc(ip,idisc) ...
                                + jacobian(this.sym.if(ip,idisc),this.strsym.xs)*M;
                        end
                    end
                end
            end
            
        case 'bdeltadisc'
            syms polydirac
            this.sym.bdeltadisc = sym(zeros(nx,ndisc));
            this.sym.bf = sym(zeros(nx,ndisc));
            this.sym.bv = sym(zeros(nx,ndisc));
            this.sym.bw = sym(zeros(nx,ndisc));
            this.sym.bz = sym(zeros(nx,ndisc));
            
            % convert deltas in xBdot
            if(ndisc>0)
                for ix = 1:nx
                    symchar = char(this.sym.xBdot(ix));
                    if(strfind(symchar,'dirac'))
                        for idisc = 1:ndisc
                            discchar = char(this.sym.rdisc(idisc));
                            str_arg_d = ['dirac(' discchar ')' ];
                            str_arg_d1 = ['dirac(1, ' discchar ')' ];
                            str_arg_1d = ['dirac(' discchar ', 1)' ];
                            str_arg_d2 = ['dirac(2, ' discchar ')' ];
                            str_arg_2d = ['dirac(' discchar ', 2)' ];
                            % find the corresponding discontinuity
                            % xdot_i = f_i(x,t) + h_i*heaviside(gdisc) + v_i*delta(gdisc) +
                            % w_i*delta(1)(gdisc) + z_i*delta(2)(gdisc);
                            % v
                            tmpstr = char(this.sym.xBdot(ix));
                            tmpstr = strrep(tmpstr,str_arg_d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.bv(ix,idisc) = cfp(2)*1/abs(this.sym.gdisc(idisc));
                            else
                                this.sym.bv(ix,idisc) = 0;
                            end
                            % w
                            tmpstr = char(this.sym.xBdot(ix));
                            tmpstr = strrep(tmpstr,str_arg_d1,'polydirac');
                            tmpstr = strrep(tmpstr,str_arg_1d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.bw(ix,idisc) = cfp(2)*this.sym.gdisc(idisc)/abs(this.sym.gdisc(idisc));
                            else
                                this.sym.bw(ix,idisc) = 0;
                            end
                            % z
                            tmpstr = char(this.sym.xBdot(ix));
                            tmpstr = strrep(tmpstr,str_arg_d2,'polydirac');
                            tmpstr = strrep(tmpstr,str_arg_2d,'polydirac');
                            cfp = coeffs(sym(tmpstr),polydirac);
                            if(length(cfp)>1)
                                this.sym.bz(ix,idisc) = cfp(2)*this.sym.gdisc(idisc)^2/abs(this.sym.gdisc(idisc));
                            else
                                this.sym.bz(ix,idisc) = 0;
                            end
                            % f
                            symchar = strrep(symchar,str_arg_d,'polydirac');
                            symchar = strrep(symchar,str_arg_d1,'polydirac');
                            symchar = strrep(symchar,str_arg_1d,'polydirac');
                            symchar = strrep(symchar,str_arg_d2,'polydirac');
                            symchar = strrep(symchar,str_arg_2d,'polydirac');
                        end
                        cfp = coeffs(sym(symchar),polydirac);
                        this.sym.bf(ix) = cfp(1);
                    else
                        this.sym.bf(ix,:) = repmat(this.sym.xBdot(ix),[1,ndisc]);
                    end
                end
                % compute bdeltadisc
                for idisc = 1:ndisc;
                    nonzero_z = any(this.sym.z(:,idisc)~=0);
                    if(nonzero_z)
                        M = jacobian(this.sym.f(:,idisc),this.strsym.xs)*this.sym.z(:,idisc) ...
                            + this.sym.w(:,idisc) ...
                            - diff(this.sym.z(:,idisc),'t') - jacobian(this.sym.z(:,idisc),this.strsym.xs)*this.sym.f(:,idisc);
                    else
                        M = this.sym.w(:,idisc);
                    end
                    nonzero_bz = any(this.sym.bz(:,idisc)~=0);
                    if(nonzero_z)
                        if(nonzero_bz)
                            N = jacobian(this.sym.bf(:,idisc),this.strsym.xs)*this.sym.z(:,idisc) ...
                                + jacobian(this.sym.bf(:,idisc),this.strsym.xBs(:))*this.sym.bz(:,idisc) ...
                                + this.sym.bw(:,idisc) ...
                                - diff(this.sym.bz(:,idisc),'t') - jacobian(this.sym.bz(:,idisc),this.strsym.xs)*this.sym.f(:,idisc) ...
                                - jacobian(this.sym.bz(:,idisc),this.strsym.xBs)*this.sym.xBdot;
                        else
                            N = jacobian(this.sym.bf(:,idisc),this.strsym.xs)*this.sym.z(:,idisc) ...
                                + this.sym.bw(:,idisc);
                        end
                    else
                        if(nonzero_bz)
                            N = jacobian(this.sym.bf(:,idisc),this.strsym.xBs(:))*this.sym.bz(:,idisc) ...
                                + this.sym.bw(:,idisc) ...
                                - diff(this.sym.bz(:,idisc),'t') - jacobian(this.sym.bz(:,idisc),this.strsym.xs)*this.sym.f(:,idisc) ...
                                - jacobian(this.sym.bz(:,idisc),this.strsym.xBs)*this.sym.xBdot;
                        else
                            N = this.sym.bw(:,idisc);
                        end
                    end
                    for ix = 1:nx
                        nonzero_M = any(M(ix,:)~=0);
                        nonzero_N = any(N(ix,:)~=0);
                        % delta_x_i =
                        % sum_j(df_idx_j*(sum_k(df_idx_k*z_k) +
                        % w_j - dot{z}_j) + v_i - dot{w}_i +
                        % ddot{z}_i
                        if(this.sym.bz(ix,idisc)~=0)
                            dbzdx = jacobian(this.sym.bz(ix,idisc),this.strsym.xs);
                            dbzdsxdot = jacobian(this.sym.bz(ix,idisc),this.strsym.xBs);
                            
                            this.sym.bdeltadisc(ix,idisc) = this.sym.bdeltadisc(ix,idisc) ...
                                + diff(this.sym.bz(ix,idisc),'t',2) ...
                                + diff(dbzdx*this.sym.xdot,'t') ...
                                + diff(dbzdsxdot*this.sym.xBdot,'t') ...
                                + jacobian(dbzdx*this.sym.f(:,idisc),this.strsym.xs)*this.sym.f(:,idisc) ...
                                + jacobian(dbzdx*this.sym.f(:,idisc),this.strsym.xBs)*this.sym.xBdot ...
                                + jacobian(dbzdsxdot*this.sym.xBdot,this.strsym.xs)*this.sym.f(:,idisc) ...
                                + jacobian(dbzdsxdot*this.sym.xBdot,this.strsym.xBs)*this.sym.xBdot;
                        end
                        
                        if(this.sym.bv(ix,idisc)~=0)
                            this.sym.bdeltadisc(ix,idisc) = this.sym.bdeltadisc(ix,idisc) ...
                                + this.sym.bv(ix,idisc);
                        end
                        
                        if(this.sym.bw(ix,idisc)~=0)
                            this.sym.bdeltadisc(ix,idisc) = this.sym.bdeltadisc(ix,idisc) ...
                                - diff(this.sym.bw(ix,idisc),'t') - jacobian(this.sym.bw(ix,idisc),this.strsym.xs)*this.sym.f(:,idisc) ...
                                - jacobian(this.sym.bw(ix,idisc),this.strsym.xBs)*this.sym.xBdot;
                            
                        end
                        
                        if(nonzero_M)
                            this.sym.bdeltadisc(ix,idisc) = this.sym.bdeltadisc(ix,idisc) ...
                                + jacobian(this.sym.bf(ix,idisc),this.strsym.xs)*M;
                        end
                        if(nonzero_N)
                            this.sym.bdeltadisc(ix,idisc) = this.sym.bdeltadisc(ix,idisc) ...
                                + jacobian(this.sym.bf(ix,idisc),this.strsym.xBs)*N;
                        end
                    end
                end
            end
            
        case 'sdeltadisc'
            
            if(ndisc>0)
                for idisc = 1:ndisc;
                    % drdp  = pdrdisc/pdp + pdrdisc/pdx*pdx/pdp
                    drdp = jacobian(this.sym.rdisc(idisc),this.strsym.ps) + jacobian(this.sym.rdisc(idisc),this.strsym.xs)*this.strsym.sx;
                    % drdt  = pdrdisc/pdt + pdrdisc/pdx*pdx/pdt
                    drdt = diff(this.sym.rdisc(idisc),'t') + jacobian(this.sym.rdisc(idisc),this.strsym.xs)*this.sym.f(:,idisc);
                    
                    % dtdp  = (1/drdt)*drdp
                    dtdp = -1/drdt*drdp;
                    
                    % dxdp  = dx/dt*dt/dp + dx/dp
                    dxdp = this.sym.f(:,idisc)*dtdp + this.strsym.sx;
                    
                    % ddelta/dp = pddelta/pdp
                    ddeltadiscdp = jacobian(this.sym.deltadisc(:,idisc),this.strsym.ps);
                    % ddelta/dt = pddelta/pdt
                    ddeltadiscdt = diff(this.sym.deltadisc(:,idisc),'t');
                    % ddelta/dx = pddelta/pdx
                    ddeltadiscdx = jacobian(this.sym.deltadisc(:,idisc),this.strsym.xs);
                    %
                    deltaxdot = mysubs(this.sym.f(:,idisc),this.strsym.xs,this.strsym.xs+this.sym.deltadisc(:,idisc))- this.sym.f(:,idisc);
                    
                    this.sym.sdeltadisc(:,:,idisc)= ...
                        - deltaxdot*dtdp ...
                        + ddeltadiscdx*dxdp ...
                        + ddeltadiscdt*dtdp ...
                        + ddeltadiscdp;
                end
            end
            
        case 'root'
            
            if(ndisc>0)
                if(isempty(this.sym.root));
                    this.sym.root = this.sym.rdisc;
                else
                    this.sym.root = [this.sym.rfun(:);this.sym.rdisc(:)];
                end
            end
            
        case 'rfun'
            this.sym.rfun = mysubs(this.sym.rfun,this.sym.p,this.strsym.ps);
            this.sym.rfun = mysubs(this.sym.rfun,this.sym.k,this.strsym.ks);
            this.sym.rfun = mysubs(this.sym.rfun,this.sym.x,this.strsym.xs);
            this.strsym.rfuns = 1;
            
        otherwise
            error('unknown function name')
    end
    
    function this = unifySyms(this,model)
        % unifySyms replaces model specific variable names with general
        % ones which conform with C naming schemes
        this.sym = mysubs(this.sym,model.sym.x,model.funs.x.sym);
        this.sym = mysubs(this.sym,model.sym.p,model.funs.p.sym);
        this.sym = mysubs(this.sym,model.sym.k,model.funs.k.sym);
    end

end