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
            this.sym = subs(this.sym,sym('D([1], am_max)'),sym('D1max'));
            this.sym = subs(this.sym,sym('D([2], am_max)'),sym('D2max'));
            this.sym = subs(this.sym,sym('am_max'),sym('max'));
        end
        
        % If we have spline, we need to parse them to get derivatives
        if (model.splineflag)
            symstr = char(this.sym);
            if (strfind(symstr, 'spline'))
                tokens = regexp(symstr, 't\,\s(\w+\.\w+)\,', 'tokens');
                nNodes = round(str2double(tokens{1}));
            end
            if (regexp(symstr, 'D\(\[(\w+|\w+\,\w+)\]\,.spline'))
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
                                this.sym(:,iCol) = subs(this.sym(:,iCol),sym(['D([' num2str(iNode*2+2) ', ' num2str(jNode*2+2) '], spline_pos)']),sym(['D' num2str(iNode*2+2) 'D' num2str(jNode*2+2) 'spline_pos']));
                                this.sym(:,iCol) = subs(this.sym(:,iCol),sym(['D([' num2str(iNode*2+2) ', ' num2str(jNode*2+2) '], spline)']),sym(['D' num2str(iNode*2+2) 'D' num2str(jNode*2+2) 'spline']));
                            end
                        end
                        this.sym(:,iCol) = subs(this.sym(:,iCol),sym(['D([' num2str(iNode*2+2) '], spline_pos)']),sym(['D' num2str(iNode*2+2) 'spline_pos']));
                        this.sym(:,iCol) = subs(this.sym(:,iCol),sym(['D([' num2str(iNode*2+2) '], spline)']),sym(['D' num2str(iNode*2+2) 'spline']));
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
        
        cstr = strrep(cstr,'log','amici::log');
        % fix derivatives again (we cant do this before as this would yield
        % incorrect symbolic expressions
        cstr = regexprep(regexprep(cstr,'D([0-9]*)([\w]*)\(','D$2\($1,'),'DD([0-9]*)([\w]*)\(','DD$2\($1,');
        cstr = strrep(strrep(cstr, 'DDspline', 'DDspline'), 'Dspline', 'Dspline');
        
        if (model.splineflag)
            if (strfind(symstr, 'spline'))
                % The floating numbers after 't' must be converted to integers
                cstr = regexprep(cstr, '([D]*(spline|spline_pos))\(t\,\w+\.\w+\,', ['amici::$1\(t\,', num2str(nNodes), '\,']);
                cstr = regexprep(cstr, '([D]*(spline|spline_pos))\((\w+)\,t\,\w+\.\w+\,', ['amici::$1\($2\,t\,', num2str(nNodes), '\,']);
                cstr = regexprep(cstr, '([D]*(spline|spline_pos))\((\w+)\,(\w+)\,t\,\w+\.\w+\,', ['amici::$1\($2\,$3\,t\,', num2str(nNodes), '\,']);
            end
        end
        
        if(numel(cstr)>1)
            
            % fix various function specific variable names/indexes
            
            cstr = regexprep(cstr,'var_x_([0-9]+)','x[$1]');
            cstr = regexprep(cstr,'var_dx_([0-9]+)','dx[$1]');
            cstr = regexprep(cstr,'var_sx_([0-9]+)','sx[$1]');
            cstr = regexprep(cstr,'var_sdx_([0-9]+)','sdx[$1]');
            % sxdot needs to preces xdot
            cstr = regexprep(cstr,'var_sxdot_([0-9]+)','sxdot[$1]');
            cstr = regexprep(cstr,'var_xdot_([0-9]+)','xdot[$1]');
            cstr = regexprep(cstr,'var_xBdot_([0-9]+)','xBdot[$1]');
            cstr = regexprep(cstr,'xdot_old_([0-9]+)','xdot_old[$1]');
            cstr = regexprep(cstr,'xdot_([0-9]+)','xdot[$1]');
            cstr = regexprep(cstr,'var_xB_([0-9]+)','xB[$1]');
            cstr = regexprep(cstr,'var_dxB_([0-9]+)','dxB[$1]');
            cstr = regexprep(cstr,'var_qBdot_([0-9]+)','qBdot[$1]');
            cstr = regexprep(cstr,'var_v_([0-9]+)', 'v[$1]');
            cstr = regexprep(cstr,'var_vB_([0-9]+)', 'vB[$1]');
            cstr = regexprep(cstr,'var_JSparse_([0-9]+)', 'JSparse->data[$1]');
            cstr = regexprep(cstr,'var_x0_([0-9]+)','x0[$1]');
            cstr = regexprep(cstr,'var_dx0_([0-9]+)','dx0[$1]');
            cstr = regexprep(cstr,'var_sx0_([0-9]+)','sx0[$1]');
            cstr = regexprep(cstr,'var_sdx0_([0-9]+)','sdx0[$1]');
            cstr = regexprep(cstr,'var_root_([0-9]+)', 'root[$1]');
            
            cstr = regexprep(cstr,'var_p_([0-9]+)','p[$1]');
            cstr = regexprep(cstr,'var_k_([0-9]+)','k[$1]');
            cstr = regexprep(cstr,'h_([0-9]+)','h[$1]');
            cstr = regexprep(cstr,'var_w_([0-9]+)','w[$1]');
            cstr = regexprep(cstr,'var_dxdotdp_([0-9]+)','dxdotdp[$1]');
            cstr = regexprep(cstr,'var_stau_([0-9]+)','stau[0]');
            cstr = regexprep(cstr,'dfdx_([0-9]+)','dfdx[$1]');
            cstr = regexprep(cstr,'M_([0-9]+)','M[$1]');
            cstr = regexprep(cstr,'var_dwdx_([0-9]+)','dwdx[$1]');
            cstr = regexprep(cstr,'var_dwdp_([0-9]+)','dwdp[$1]');
            cstr = regexprep(cstr,'tmp_J_([0-9]+)','J->data[$1]');
            cstr = regexprep(cstr,'tmp_dxdotdp_([0-9]+)','dxdotdp[$1]');
            
            cstr = regexprep(cstr,'var_y_([0-9]+)','y[$1]');
            cstr = regexprep(cstr,'my_([0-9]+)','my[$1]');
            cstr = regexprep(cstr,'var_z_([0-9]+)','z[$1]');
            cstr = regexprep(cstr,'mz_([0-9]+)','mz[$1]');
            cstr = regexprep(cstr,'var_rz_([0-9]+)','rz[$1]');
            cstr = regexprep(cstr,'var_srz_([0-9]+)','srz[$1]');
            cstr = regexprep(cstr,'var_sy_([0-9]+)','sy[$1]');
            cstr = regexprep(cstr,'var_sz_([0-9]+)','sz[$1]');
           
            cstr = regexprep(cstr,'var_dydx[_\[]*([0-9\+\*]+)[\]]*','dydx[$1]'); % matches both _... and [...]
            cstr = regexprep(cstr,'var_dzdx[_\[]*([0-9\+\*]+)[\]]*','dzdx[$1]');
            cstr = regexprep(cstr,'var_drzdx[_\[]*([0-9\+\*]+)[\]]*','drzdx[$1]');
            cstr = regexprep(cstr,'var_dydp_([0-9]+)',['dydp[$1]']);
            cstr = regexprep(cstr,'var_dzdp_([0-9]+)',['dzdp[$1]']);
            cstr = regexprep(cstr,'var_drzdp_([0-9]+)',['drzdp[$1]']);
            cstr = regexprep(cstr,'var_drootdp_([0-9]+)',['drzdp[$1]']);
            cstr = regexprep(cstr,'var_deltax_([0-9]+)','deltax[$1]');
            cstr = regexprep(cstr,'var_deltaxB_([0-9]+)','deltaxB[$1]');
            cstr = regexprep(cstr,'var_deltasx_([0-9]+)','deltasx[$1]');
            cstr = regexprep(cstr,'var_deltaqB_([0-9]+)','deltaqB[$1]');
            cstr = regexprep(cstr,'var_sigma_y_([0-9]+)','sigmay[$1]');
            cstr = regexprep(cstr,'var_sigma_z_([0-9]+)','sigmaz[$1]');
            cstr = regexprep(cstr,'var_dsigma_zdp_([0-9]+)',['dsigmazdp[$1]']);
            cstr = regexprep(cstr,'var_dsigma_ydp_([0-9]+)',['dsigmaydp[$1]']);
            
            cstr = regexprep(cstr,'var_dsdydp_([0-9]+)',['dsigmaydp[ip*' num2str(model.ny) ' + $1]']);
            cstr = regexprep(cstr,'var_dsdzdp_([0-9]+)',['dsigmazdp[ip*' num2str(model.nz) ' + $1]']);
            cstr = regexprep(cstr,'var_Jy_([0-9]+)','nllh[$1]');
            cstr = regexprep(cstr,'var_Jz_([0-9]+)','nllh[$1]');
            cstr = regexprep(cstr,'var_Jrz_([0-9]+)','nllh[$1]'); % not a typo; we dont want to creat an additional variable that we end the end add to Jz anyways.
            cstr = regexprep(cstr,'var_dJydy[_\[]*([0-9\+\*]+)[\]]*','dJydy[$1]'); % matches both _... and [...]
            cstr = regexprep(cstr,'var_dJzdz[_\[]*([0-9\+\*]+)[\]]*','dJzdz[$1]');
            cstr = regexprep(cstr,'var_dJrzdz[_\[]*([0-9\+\*]+)[\]]*','dJrzdz[$1]');
            cstr = regexprep(cstr,'var_dJydsigma[_\[]*([0-9\+\*]+)[\]]*','dJydsigma[$1]');
            cstr = regexprep(cstr,'var_dJzdsigma[_\[]*([0-9\+\*]+)[\]]*','dJzdsigma[$1]');
            cstr = regexprep(cstr,'var_dJrzdsigma[_\[]*([0-9\+\*]+)[\]]*','dJrzdsigma[$1]');
            cstr = regexprep(cstr,'var_JDiag[_\[]*([0-9\+\*]+)[\]]*','JDiag[$1]');
        end
        
        %%
        % print to file
        fprintf(fid,[cstr '\n']);
    end
end
