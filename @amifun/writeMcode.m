function writeMcode(this,model)
% writeMcode generates matlab evaluable code for specific model functions
%
% Parameters:
%  model: model defintion object @type amimodel
%
% Return values:
%  void


% replace heaviside placeholder vars
if(~isempty(model.event))
    h_vars = [sym('h_0');sym('h_',[length(model.event)-1,1])];
    h_rep= sym(zeros(size(h_vars)));
    model.getFun([],'root');
    for ievent = 1:length(model.event)
        h_rep(ievent) = heaviside(model.fun.root.sym(ievent));
    end
    this.sym = subs(this.sym,h_vars,h_rep);
end
    

if(strcmp(this.funstr,'w') || strcmp(this.funstr,'dwdp') || strcmp(this.funstr,'dwdx'))
    fid = fopen(fullfile(model.wrap_path,'models',model.modelname,['calc_' this.funstr '_',model.modelname,'.m']), 'w');
    fprintf(fid,['%% calc_' this.funstr  '_' model.modelname '.m is an additional matlab file to\n'...
        '%% compute helping variables for the jacobian with AMICI.\n\n']);
    fprintf(fid, ['function ' this.funstr ' = calc_' this.funstr  '_' model.modelname '(p, x, k, t)\n\n']);
    
    % Assignment of inputs
    for l = 1 : numel(model.sym.p)
        fprintf(fid, ['p_' num2str(l-1) ' = p(' num2str(l) ');\n']);
    end
    for l = 1 : numel(model.sym.x)
        fprintf(fid, ['x_' num2str(l-1) ' = x(' num2str(l) ');\n']);
    end
    for l = 1 : numel(model.sym.k)
        fprintf(fid, ['k_' num2str(l-1) ' = k(' num2str(l) ');\n']);
    end
    fprintf(fid, '\n');

    % Writing and giving output
    for l = 1 : numel(this.sym)
        fprintf(fid, [this.funstr '_' num2str(l-1) ' = ' char(this.sym(l)) ';\n']);
    end
    fprintf(fid,[this.funstr ' = [']);
    for l = 1 : numel(this.sym)
        if (l==1)
            fprintf(fid, [this.funstr '_' num2str(l-1) ]);
        else
            fprintf(fid, [', ' this.funstr '_' num2str(l-1) ]);
            if (abs((l/10.0) - round(l/10.0)) < 0.001)
                fprintf(fid, ' ...\n    ');
            end
        end
    end
    fprintf(fid, '];\n\n');
    fprintf(fid, 'end');
    fclose(fid);
else
    mfun(this.sym, 'file', fullfile(model.wrap_path,'models',...
        model.modelname,['eval_' this.funstr '_',model.modelname,'.m']), ...
        'vars', {model.fun.p.sym, model.fun.x.sym, model.fun.k.sym, ...
        model.fun.w.strsym(find(model.fun.w.strsym)), ...
        model.fun.dwdp.strsym(find(model.fun.dwdp.strsym)), ...
        model.fun.dwdx.strsym(find(model.fun.dwdx.strsym)), 't'});
end


end