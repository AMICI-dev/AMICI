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
    if(~isfield(model.fun,'root'))
        model.getFun([],'root');
    end
    for ievent = 1:length(model.event)
        h_rep(ievent) = heaviside(model.fun.root.sym(ievent));
    end
    this.sym_noopt = subs(this.sym_noopt,h_vars,h_rep);
end
    
ami_mfun(this.sym_noopt, 'file', fullfile(model.wrap_path,'models',...
    model.modelname,[ this.funstr '_',model.modelname,'.m']), ...
    'vars', {'t',model.fun.x.sym,model.fun.p.sym,model.fun.k.sym},'varnames',{'t','x','p','k'});
end