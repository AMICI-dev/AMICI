function generateJacobianFile(this)
    % generateJacobianFile generates files to compute the Jacobian.

    % Create the files tp create w, dwdwp, dwdx
    fid = fopen(fullfile(this.wrap_path,'models',this.modelname,['calc_w_',this.modelname,'.m']), 'w');
    fprintf(fid,['%% calc_w_' this.modelname '.m is an additional matlab file to\n'...
        '%% compute helping variables for the jacobian with AMICI.\n\n']);
    fprintf(fid, ['function w = calc_w_' this.modelname '(p, x, k, t)\n\n']);

    % Assignment of inputs
    for l = 1 : numel(this.sym.p)
        fprintf(fid, ['p_' num2str(l-1) ' = p(' num2str(l) ');\n']);
    end
    for l = 1 : numel(this.sym.x)
        fprintf(fid, ['x_' num2str(l-1) ' = x(' num2str(l) ');\n']);
    end
    for l = 1 : numel(this.sym.k)
        fprintf(fid, ['k_' num2str(l-1) ' = k(' num2str(l) ');\n']);
    end
    fprintf(fid, '\n');

    % Writing and giving output
    for l = 1 : numel(this.fun.w.sym)
        fprintf(fid, ['w_' num2str(l-1) ' = ' char(this.fun.w.sym(l)) ';\n']);
    end
    fprintf(fid, 'w = [');
    for l = 1 : numel(this.fun.w.sym)
        if (l==1)
            fprintf(fid, ['w_' num2str(l-1) ]);
        else
            fprintf(fid, [', w_' num2str(l-1) ]);
            if (abs((l/10.0) - round(l/10.0)) < 0.001)
                fprintf(fid, ' ...\n    ');
            end
        end
    end
    fprintf(fid, '];\n\n');
    fprintf(fid, 'end');
    fclose(fid);




    % Writing the dwdp file
    fid = fopen(fullfile(this.wrap_path,'models',this.modelname,['calc_dwdp_',this.modelname,'.m']),'w');
    fprintf(fid,['%% calc_dwdp_' this.modelname '.m is an additional matlab file to\n'...
        '%% compute helping variables for the jacobian with AMICI.\n\n']);
    fprintf(fid, ['function dwdp = calc_dwdp_' this.modelname '(p, x, k, w, t)\n\n']);

    % Assignment of inputs
    for l = 1 : numel(this.sym.p)
        fprintf(fid, ['p_' num2str(l-1) ' = p(' num2str(l) ');\n']);
    end
    for l = 1 : numel(this.sym.x)
        fprintf(fid, ['x_' num2str(l-1) ' = x(' num2str(l) ');\n']);
    end
    for l = 1 : numel(this.sym.k)
        fprintf(fid, ['k_' num2str(l-1) ' = k(' num2str(l) ');\n']);
    end
    for l = 1 : numel(this.fun.w.sym)
        fprintf(fid, ['w_' num2str(l-1) ' = w(' num2str(l) ');\n']);
    end
    fprintf(fid, '\n');


    % Writing and giving output
    for l = 1 : numel(this.fun.dwdp.sym)
        fprintf(fid, ['dwdp_' num2str(l-1) ' = ' char(this.fun.dwdp.sym(l)) ';\n']);
    end
    fprintf(fid, 'dwdp = [');
    for l = 1 : numel(this.fun.dwdp.sym)
        if (l==1)
            fprintf(fid, ['dwdp_' num2str(l-1) ]);
        else
            fprintf(fid, [', dwdp_' num2str(l-1) ]);
            if (abs((l/10.0) - round(l/10.0)) < 0.001)
                fprintf(fid, ' ...\n    ');
            end
        end
    end
    fprintf(fid, '];\n\n');
    fprintf(fid, 'end');
    fclose(fid);


    % Writing the dwdx file
    fid = fopen(fullfile(this.wrap_path,'models',this.modelname,['calc_dwdx_',this.modelname,'.m']),'w');
    fprintf(fid,['%% calc_dwdx_' this.modelname '.m is an additional matlab file to\n'...
        '%% compute helping variables for the jacobian with AMICI.\n\n']);
    fprintf(fid, ['function dwdx = calc_dwdx_' this.modelname '(p, x, k, w, t)\n\n']);

    % Assignment of inputs
    for l = 1 : numel(this.sym.p)
        fprintf(fid, ['p_' num2str(l-1) ' = p(' num2str(l) ');\n']);
    end
    for l = 1 : numel(this.sym.x)
        fprintf(fid, ['x_' num2str(l-1) ' = x(' num2str(l) ');\n']);
    end
    for l = 1 : numel(this.sym.k)
        fprintf(fid, ['k_' num2str(l-1) ' = k(' num2str(l) ');\n']);
    end
    for l = 1 : numel(this.fun.w.sym)
        fprintf(fid, ['w_' num2str(l-1) ' = w(' num2str(l) ');\n']);
    end
    fprintf(fid, '\n');

    % Writing and giving output
    for l = 1 : numel(this.fun.dwdx.sym)
        fprintf(fid, ['dwdx_' num2str(l-1) ' = ' char(this.fun.dwdx.sym(l)) ';\n']);
    end
    fprintf(fid, 'dwdx = [');
    for l = 1 : numel(this.fun.dwdx.sym)
        if (l==1)
            fprintf(fid, ['dwdx_' num2str(l-1) ]);
        else
            fprintf(fid, [', dwdx_' num2str(l-1) ]);
            if (abs((l/10.0) - round(l/10.0)) < 0.001)
                fprintf(fid, ' ...\n    ');
            end
        end
    end
    fprintf(fid, '];\n\n');
    fprintf(fid, 'end');
    fclose(fid);


    mfun(this.fun.J.sym, 'file', fullfile(this.wrap_path,'models',...
        this.modelname,['eval_J_',this.modelname,'.m']), ...
        'vars', {this.fun.p.sym, this.fun.x.sym, this.fun.k.sym, ...
        this.fun.w.strsym(this.fun.w.strsym~=0), ...
        this.fun.dwdp.strsym(this.fun.dwdp.strsym~=0), ...
        this.fun.dwdx.strsym(this.fun.dwdx.strsym~=0), 't'});
    mfun(this.fun.dxdotdp.sym, 'file', fullfile(this.wrap_path, ...
        'models',this.modelname,['eval_dxdotdp_',this.modelname,'.m']), ...
        'vars', {this.fun.p.sym, this.fun.x.sym, this.fun.k.sym, ...
        this.fun.w.strsym(this.fun.w.strsym~=0), ...
        this.fun.dwdp.strsym(this.fun.dwdp.strsym~=0), ...
        this.fun.dwdx.strsym(this.fun.dwdx.strsym~=0), 't'});    
end

