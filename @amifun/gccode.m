function this = gccode(this,model,fid)
    % gccode transforms symbolic expressions into c code and writes the
    % respective expression into a specified file
    %
    % Parameters:
    %  model: model defintion object @type amimodel
    %  fid: file id in which the expression should be written @type fileid
    %
    % Return values:
    %  this: function definition object @type amifun
    
    
    if(any(any(any(this.sym~=0))))
        
        % replace unknown partial derivatives
        if(model.maxflag)
            this.sym = subs(this.sym,sym('D([1], am_max)'),sym('D1am_max'));
            this.sym = subs(this.sym,sym('D([2], am_max)'),sym('D2am_max'));
        end
        if(model.minflag)
            this.sym = subs(this.sym,sym('D([1], am_min)'),sym('D1am_min'));
            this.sym = subs(this.sym,sym('D([2], am_min)'),sym('D2am_min'));
        end

        if(model.splineflag)
            for nodes = [3,4,5,10]
                for ideriv = 1:nodes
                    this.sym = subs(this.sym,sym(['D([' num2str(ideriv*2+1) '], spline_pos' num2str(nodes) ')']),sym(['D' num2str(ideriv*2+1) 'spline_pos' num2str(nodes)]));
                    this.sym = subs(this.sym,sym(['D([' num2str(ideriv*2+1) '], spline' num2str(nodes) ')']),sym(['D' num2str(ideriv*2+1) 'spline' num2str(nodes)]));
                end
            end
        end
        
        cstr = ccode(this.sym);
        if(numel(this.sym)>1)
            if(strcmp(this.funstr,'J') || strcmp(this.funstr,'JB'))
                cstr = regexprep(cstr,'T\[([0-9]*)\]\[([0-9]*)\]',[this.cvar '[$1*' num2str(model.nx) '+$2]']);
            else
                cstr = strrep(cstr,'T[0]',this.cvar);
            end
        else
            cstr = strrep(cstr,'t0',[this.cvar '[0]']);
        end
        cstr = strrep(cstr,'log','amilog');
        % fix derivatives again (we cant do this before as this would yield
        % incorrect symbolic expressions
        cstr = regexprep(cstr,'D([0-9]*)([\w]*)\(','D$2\($1,');
        
        
        if(numel(cstr)>1)
            
            % fix various function specific variable names/indexes
            
            if(ismember('x',this.nvecs'))
                cstr = regexprep(cstr,'x_([0-9]+)','x\[$1\]');
            else
                if(ismember('*sx',this.nvecs'))
                    cstr = regexprep(cstr,'sx_([0-9]+)','sx\[$1\]');
                end
            end
            if(ismember('xdot',this.nvecs'))
                cstr = regexprep(cstr,'xdot_([0-9]+)','xdot\[$1\]');
            end
            if(ismember('xdot_old',this.nvecs'))
                cstr = regexprep(cstr,'xdot_old_([0-9]+)','xdot_old\[$1\]');
            end
            cstr = regexprep(cstr,'dwdp_([0-9]+)','dwdp_tmp\[$1\]');
            cstr = regexprep(cstr,'dydp_([0-9]+)','dydp\[$1\]');
            cstr = regexprep(cstr,'p_([0-9]+)','p\[$1\]');
            cstr = regexprep(cstr,'k_([0-9]+)','k\[$1\]');
            cstr = regexprep(cstr,'h_([0-9]+)','h\[$1\]');
            cstr = regexprep(cstr,'w_([0-9]+)','w_tmp\[$1\]');
            
            cstr = regexprep(cstr,'dwdx_([0-9]+)','dwdx_tmp\[$1\]');
            cstr = regexprep(cstr,'dydx_([0-9]+)','dydx\[$1\]');
            if(ismember('xB',this.nvecs'))
                cstr = regexprep(cstr,'xB_([0-9]+)','xB\[$1\]');
            end
            cstr = regexprep(cstr,'J([0-9]+)','J\[$1\]');
            cstr = regexprep(cstr,'xdotdp([0-9]+)',['xdotdp\[$1 + ip*' num2str(model.nx) '\]']);
            
            if(strcmp(this.cvar,'qBdot'))
                cstr = regexprep(cstr,'qBdot\[([0-9]*)\]','qBdot\[ip]');
            elseif(strcmp(this.cvar,'y'))
                cstr = regexprep(cstr,[this.cvar '\[([0-9]*)\]'],[this.cvar '\[it+nt*$1\]']);
            elseif(strcmp(this.cvar,'sy'))
                cstr = regexprep(cstr,[this.cvar '\[([0-9]*)\]'],[this.cvar '\[it+nt*\($1+ip*' num2str(model.ny) '\)\]']);
            elseif(strcmp(this.cvar,'z'))
                cstr = regexprep(cstr,[this.cvar '\[([0-9]*)\]'],[this.cvar '\[nroots[ie] + nmaxevent*$1\]']);
            elseif(strcmp(this.cvar,'sz'))
                cstr = regexprep(cstr,[this.cvar '\[([0-9]*)\]'],[this.cvar '\[nroots[ie] + nmaxevent*($1+ip*' num2str(model.nz) ')\]']);
            elseif(strcmp(this.cvar,'dxdotdp'))
                cstr = regexprep(cstr, 'dxdotdp\[([0-9]*)\]',['dxdotdp[\($1+ip*' num2str(model.nx) '\)\]']);
            elseif(strcmp(this.cvar,'Jdata'))
                cstr = strrep(cstr, 'Jdata[', 'J->data[');
            elseif(strcmp(this.cvar,'JBdata'))
                cstr = strrep(cstr, 'JBdata[', 'JB->data[');
            elseif(strcmp(this.cvar,'Jv'))
                cstr = regexprep(cstr, 'v_([0-9]*)', 'v_tmp\[$1\]');
            elseif(strcmp(this.cvar,'JvB'))
                cstr = regexprep(cstr, 'v_([0-9]*)', 'vB_tmp\[$1\]');
            elseif(strcmp(this.cvar,'sxdot'))
                cstr = strrep(cstr,'tmp_J[','tmp_J->data[');
            elseif(strcmp(this.cvar,'dsigma_zdp')  || strcmp(this.cvar,'sz')|| strcmp(this.cvar,'dzdp'))
                cstr = regexprep(cstr, [this.cvar '\[([0-9]*)\] = '], [ this.cvar '[ip*' num2str(model.nz) ' + $1] = ']);
            elseif(strcmp(this.cvar,'dsigma_ydp') || strcmp(this.cvar,'dydp'))
                cstr = regexprep(cstr, [this.cvar '\[([0-9]*)\] = '], [ this.cvar '[ip*' num2str(model.ny) ' + $1] = ']);
            elseif(strcmp(this.cvar,'deltasx'))
                cstr = regexprep(cstr, [this.cvar '\[([0-9]*)\] = '], [ this.cvar '[ip*' num2str(model.nx) ' + $1] = ']);
            elseif(strcmp(this.cvar,'dJydx'))
                cstr = regexprep(cstr, [this.cvar '\[([0-9]*)\] = '], [ this.cvar '[it+$1*nt] = ']);
            elseif(strcmp(this.cvar,'dJzdx'))
                cstr = regexprep(cstr, [this.cvar '\[([0-9]*)\] = '], [ this.cvar '[nroots[ie]+$1*maxevent] = ']);
            end
            
            if(strfind(this.cvar,'Jy'))
                cstr = regexprep(cstr,'dydx_([0-9]+)','dydx\[$1]');
                cstr = regexprep(cstr,'dydp_([0-9]+)',['dydp\[$1+ip*' num2str(model.ny) ']']);
                cstr = regexprep(cstr,'my_([0-9]+)','my\[it+nt*$1]');            
                cstr = regexprep(cstr,'sdy_([0-9]+)','sd_y\[$1\]');
                cstr = regexprep(cstr,'dsdydp\[([0-9]*)\]','dsigma_ydp\[$1\]');
                if(strcmp(this.cvar,'sJy'))
                    cstr = regexprep(cstr,'sy_([0-9]+)','sy\[$1\]');
                    cstr = regexprep(cstr,'dJydy_([0-9]+)','dJydy\[$1\]');
                else
                    cstr = regexprep(cstr,'sy_([0-9]+)',['sy\[it+nt*\($1+ip*' num2str(model.ny) '\)\]']);
                    cstr = regexprep(cstr,'y_([0-9]+)','y\[it+nt*$1\]');
                    cstr = strrep(cstr,'=','+=');
                end
            end
            
            if(strfind(this.cvar,'Jz'))
                cstr = regexprep(cstr,'dzdx_([0-9]+)','dzdx\[$1]');
                cstr = regexprep(cstr,'dzdp_([0-9]+)',['dzdp\[$1+ip*' num2str(model.nz) ']']);
                cstr = regexprep(cstr,'mz_([0-9]+)','mz\[nroots[ie]+nmaxevent*$1]');
                cstr = regexprep(cstr,'sz_([0-9]+)',['sz\[nroots[ie]+nmaxevent*\($1+ip*' num2str(model.nz) '\)\]']);
                cstr = regexprep(cstr,'sdz_([0-9]+)','sd_z\[$1\]');
                cstr = regexprep(cstr,'dsdzdp\[([0-9]*)\]','dsigma_zdp\[$1\]');
                cstr = regexprep(cstr,'z_([0-9]+)','z\[nroots[ie]+nmaxevent*$1\]');
                cstr = strrep(cstr,'=','+=');
            end
            
            
            
            for nvec = this.nvecs
                nstr = strrep(nvec{1},'*','');
                cstr = strrep(cstr, [nstr '['],[nstr '_tmp[']);
            end
            if(strcmp(this.funstr,'dydx'))
                cstr = strrep(cstr,'dydx_tmp','dydx');
            end
            if(strcmp(this.funstr,'deltax'))
                cstr = strrep(cstr,'deltax_tmp','deltax');
            end
            if(strcmp(this.funstr,'deltasx'))
                cstr = strrep(cstr,'deltasx_tmp','deltasx');
            end
            if(strcmp(this.funstr,'deltaxB'))
                cstr = strrep(cstr,'deltaxB_tmp','deltaxB');
            end
            if(strcmp(this.funstr,'deltaqB'))
                cstr = strrep(cstr,'deltaqB_tmp','deltaqB');
            end
            if(strcmp(this.funstr,'dJydx'))
                cstr = strrep(cstr,'ydx_tmp','ydx');
            end
            
            
            %%
            % print to file
            fprintf(fid,[cstr '\n']);
        end
    end
end