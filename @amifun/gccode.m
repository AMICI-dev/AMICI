function this = gccode(this,model,fid)
    % gccode transforms symbolic expressions into c code and writes the
    % respective expression into a specified file
    %
    % Parameters:
    %  model: model definition object @type amimodel
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
        
        % If we have spline, we need to parse them to get derivatives
        if (model.splineflag)
            symstr = char(this.sym);
            if (strfind(symstr, 'spline'))
                tokens = regexp(symstr, 't\,\s(\w+\.\w+)\,', 'tokens');
                nNodes = round(str2double(tokens{1}));
            end
            if (regexp(symstr, 'D\(\[(\w+|\w+\,\w+)\]\,.am_spline'))
                isDSpline = true;
            else
                isDSpline = false;
            end
            
            if (isDSpline)
                [~, nCol] = size(this.sym);
                for iCol = 1 : nCol
                    for iNode = 1 : nNodes
                        if (model.o2flag)
                            for jNode = 1:nNodes
                                this.sym(:,iCol) = subs(this.sym(:,iCol),sym(['D([' num2str(iNode*2+2) ', ' num2str(jNode*2+2) '], am_spline_pos)']),sym(['D' num2str(iNode*2+2) 'D' num2str(jNode*2+2) 'am_spline_pos']));
                                this.sym(:,iCol) = subs(this.sym(:,iCol),sym(['D([' num2str(iNode*2+2) ', ' num2str(jNode*2+2) '], am_spline)']),sym(['D' num2str(iNode*2+2) 'D' num2str(jNode*2+2) 'am_spline']));
                            end
                        end
                        this.sym(:,iCol) = subs(this.sym(:,iCol),sym(['D([' num2str(iNode*2+2) '], am_spline_pos)']),sym(['D' num2str(iNode*2+2) 'am_spline_pos']));
                        this.sym(:,iCol) = subs(this.sym(:,iCol),sym(['D([' num2str(iNode*2+2) '], am_spline)']),sym(['D' num2str(iNode*2+2) 'am_spline']));
                    end
                end
            end
        end
        
        cstr = ccode(this.sym);
        if(~strcmp(cstr(3:4),'t0'))
            if(any(strcmp(this.funstr,{'J','JB','JDiag','dJydsigma','dJydy','dJzdsigma','dJzdz','dJrzdsigma','dJrzdz','dydx','dzdx','drzdx','M','dfdx'}) ))
                cstr = regexprep(cstr,'T\[([0-9]*)\]\[([0-9]*)\]',[this.cvar '[$1+$2*' num2str(size(this.sym,1)) ']']);
            else
                cstr = regexprep(cstr,'T\[([0-9]*)\]\[0\]',[this.cvar '_$1']);
                cstr = regexprep(cstr,'T\[0\]\[([0-9]*)\]',[this.cvar '_$1']);
            end
        else
            cstr = strrep(cstr,'t0',[this.cvar '_0']);
        end
        
        cstr = strrep(cstr,'log','amilog');
        % fix derivatives again (we cant do this before as this would yield
        % incorrect symbolic expressions
        cstr = regexprep(regexprep(cstr,'D([0-9]*)([\w]*)\(','D$2\($1,'),'DD([0-9]*)([\w]*)\(','DD$2\($1,');
        cstr = strrep(strrep(cstr, 'DDam_spline', 'am_DDspline'), 'Dam_spline', 'am_Dspline');
        
        if (model.splineflag)
            if (strfind(symstr, 'spline'))
                % The floating numbers after 't' must be converted to integers
                cstr = regexprep(cstr, '(spline|spline_pos)\(t\,\w+\.\w+\,', ['$1\(t\,', num2str(nNodes), '\,']);
                cstr = regexprep(cstr, '(spline|spline_pos)\((\w+)\,t\,\w+\.\w+\,', ['$1\($2\,t\,', num2str(nNodes), '\,']);
                cstr = regexprep(cstr, '(spline|spline_pos)\((\w+)\,(\w+)\,t\,\w+\.\w+\,', ['$1\($2\,$3\,t\,', num2str(nNodes), '\,']);
            end
        end
        
        if(numel(cstr)>1)
            
            % fix various function specific variable names/indexes
            
            cstr = regexprep(cstr,'var_x_([0-9]+)','x_tmp[$1]');
            cstr = regexprep(cstr,'var_sx_([0-9]+)','sx_tmp[$1]');
            % sxdot needs to preces xdot
            cstr = regexprep(cstr,'var_sxdot_([0-9]+)','sxdot_tmp[$1]');
            cstr = regexprep(cstr,'var_xdot_([0-9]+)','xdot_tmp[$1]');
            cstr = regexprep(cstr,'var_xBdot_([0-9]+)','xBdot_tmp[$1]');
            cstr = regexprep(cstr,'xdot_old_([0-9]+)','xdot_old_tmp[$1]');
            cstr = regexprep(cstr,'xdot_([0-9]+)','xdot_tmp[$1]');
            cstr = regexprep(cstr,'var_xB_([0-9]+)','xB_tmp[$1]');
            cstr = regexprep(cstr,'var_qBdot_([0-9]+)','qBdot_tmp[ip + udata->nplist*$1]');
            cstr = regexprep(cstr,'var_v_([0-9]+)', 'v_tmp[$1]');
            cstr = regexprep(cstr,'var_vB_([0-9]+)', 'vB_tmp[$1]');
            cstr = strrep(cstr, 'Jdata[', 'J->data[');
            cstr = strrep(cstr, 'JBdata[', 'JB->data[');
            cstr = regexprep(cstr, 'Jdata_([0-9]+)', 'J->data[$1]');
            cstr = regexprep(cstr, 'JBdata_([0-9]+)', 'JB->data[$1]');
            cstr = regexprep(cstr, 'Jv_tmp_([0-9]+)', 'Jv_tmp[$1]');
            cstr = regexprep(cstr, 'JvB_tmp_([0-9]+)', 'JvB_tmp[$1]');
            cstr = regexprep(cstr,'var_x0_([0-9]+)','x0_tmp[$1]');
            cstr = regexprep(cstr,'var_dx0_([0-9]+)','dx0_tmp[$1]');
            cstr = regexprep(cstr,'var_sx0_([0-9]+)','sx0_tmp[$1]');
            cstr = regexprep(cstr,'var_sdx0_([0-9]+)','sdx0_tmp[$1]');
            cstr = regexprep(cstr,'var_root_([0-9]+)', 'root[$1]');
            
            cstr = regexprep(cstr,'var_p_([0-9]+)','udata->p[$1]');
            cstr = regexprep(cstr,'var_k_([0-9]+)','udata->k[$1]');
            cstr = regexprep(cstr,'h_([0-9]+)','udata->h[$1]');
            cstr = regexprep(cstr,'var_w_([0-9]+)','udata->w[$1]');
            cstr = regexprep(cstr,'var_dxdotdp_([0-9]+)','udata->dxdotdp[$1 + ip*udata->nx]');
            cstr = regexprep(cstr,'var_stau_([0-9]+)','udata->stau[ip]');
            cstr = regexprep(cstr,'var_dwdx_([0-9]+)','udata->dwdx[$1]');
            cstr = regexprep(cstr,'var_dwdp_([0-9]+)','udata->dwdp[$1]');
            cstr = regexprep(cstr,'tmp_J_([0-9]+)','udata->J->data[$1]');
            cstr = regexprep(cstr,'tmp_dxdotdp_([0-9]+)','udata->dxdotdp[$1 + ip*udata->nx]');
            
            cstr = regexprep(cstr,'var_y_([0-9]+)','rdata->y[it + udata->nt*$1]');
            cstr = regexprep(cstr,'var_z_([0-9]+)','rdata->z[tdata->nroots[ie]+udata->nmaxevent*$1]');
            cstr = regexprep(cstr,'var_rz_([0-9]+)','rdata->rz[tdata->nroots[ie]+udata->nmaxevent*$1]');
            cstr = regexprep(cstr,'var_srz_([0-9]+)','rdata->srz[tdata->nroots[ie]+udata->nmaxevent*(ip*udata->nz + $1)]');
            cstr = regexprep(cstr,'var_sy_([0-9]+)','rdata->sy[it+ udata->nt*($1+ip * udata->ny)]');
            cstr = regexprep(cstr,'var_sz_([0-9]+)','rdata->sz[tdata->nroots[ie] + udata->nmaxevent*(ip*udata->nz + $1)]');
           
            cstr = regexprep(cstr,'var_dydx[_\[]*([0-9\+\*]+)[\]]*','tdata->dydx[$1]'); % matches both _... and [...]
            cstr = regexprep(cstr,'var_dzdx[_\[]*([0-9\+\*]+)[\]]*','tdata->dzdx[$1]');
            cstr = regexprep(cstr,'var_drzdx[_\[]*([0-9\+\*]+)[\]]*','tdata->drzdx[$1]'); 
            cstr = regexprep(cstr,'var_dydp_([0-9]+)','tdata->dydp[ip*udata->ny + $1]');    
            cstr = regexprep(cstr,'var_dzdp_([0-9]+)','tdata->dzdp[ip*udata->nz + $1]');
            cstr = regexprep(cstr,'var_drzdp_([0-9]+)','tdata->drzdp[ip*udata->nz + $1]');
            cstr = regexprep(cstr,'var_drootdp_([0-9]+)','tdata->drzdp[ip*udata->nz + $1]');
            cstr = regexprep(cstr,'var_deltax_([0-9]+)','tdata->deltax[$1]');
            cstr = regexprep(cstr,'var_deltaxB_([0-9]+)','tdata->deltaxB[$1]');
            cstr = regexprep(cstr,'var_deltasx_([0-9]+)','tdata->deltasx[ip*udata->nx + $1]');
            cstr = regexprep(cstr,'var_deltaqB_([0-9]+)','tdata->deltaqB[ip+$1]');
            cstr = regexprep(cstr,'var_sigma_y_([0-9]+)','tdata->sigmay[$1]');
            cstr = regexprep(cstr,'var_sigma_z_([0-9]+)','tdata->sigmaz[$1]');
            cstr = regexprep(cstr,'var_dsigma_zdp_([0-9]+)','tdata->dsigmazdp[ip*udata->nz + $1]');
            cstr = regexprep(cstr,'var_dsigma_ydp_([0-9]+)','tdata->dsigmaydp[ip*udata->ny + $1]');
            
            cstr = regexprep(cstr,'var_dsdydp_([0-9]+)','tdata->dsigmaydp[ip*udata->ny + $1]');
            cstr = regexprep(cstr,'var_dsdzdp_([0-9]+)','tdata->dsigmazdp[ip*udata->nz + $1]');
            cstr = regexprep(cstr,'var_Jy_([0-9]+)','tdata->Jy[$1]');
            cstr = regexprep(cstr,'var_Jz_([0-9]+)','tdata->Jz[$1]');
            cstr = regexprep(cstr,'var_Jrz_([0-9]+)','tdata->Jz[$1]'); % not a typo; we dont want to creat an additional variable that we end the end add to Jz anyways.
            cstr = regexprep(cstr,'var_dJydy[_\[]*([0-9\+\*]+)[\]]*','tdata->dJydy[iy+($1)*udata->nytrue]'); % matches both _... and [...]
            cstr = regexprep(cstr,'var_dJzdz[_\[]*([0-9\+\*]+)[\]]*','tdata->dJzdz[iz+($1)*udata->nztrue]'); 
            cstr = regexprep(cstr,'var_dJrzdz[_\[]*([0-9\+\*]+)[\]]*','tdata->dJrzdz[iz+($1)*udata->nztrue]');
            cstr = regexprep(cstr,'var_dJydsigma[_\[]*([0-9\+\*]+)[\]]*','tdata->dJydsigma[iy+($1)*udata->nytrue]');
            cstr = regexprep(cstr,'var_dJzdsigma[_\[]*([0-9\+\*]+)[\]]*','tdata->dJzdsigma[iz+($1)*udata->nztrue]');
            cstr = regexprep(cstr,'var_dJrzdsigma[_\[]*([0-9\+\*]+)[\]]*','tdata->dJrzdsigma[iz+($1)*udata->nztrue]');
           
            
            if(~isempty(strfind(this.cvar,'Jy')))
                cstr = regexprep(cstr,'my_([0-9]+)','edata->my[it+udata->nt*$1]');
            end
            if(~isempty([strfind(this.cvar,'Jz'),strfind(this.cvar,'Jrz')]))
                cstr = regexprep(cstr,'mz_([0-9]+)','edata->mz[tdata->nroots[ie]+udata->nmaxevent*$1]');
            end
            if(ismember(this.cvar,{'var_Jy','var_Jz','var_Jrz'}))
                cstr = strrep(cstr,'=','+=');
            end
        end
        
        %%
        % print to file
        fprintf(fid,[cstr '\n']);
    end
end