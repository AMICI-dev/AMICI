function compileC(this)
    % compileC compiles the mex simulation file
    %
    % Return values:
    %  this: model definition object @type amimodel

    amimodel.compileAndLinkModel(this.modelname, this.wrap_path, this.recompile, this.coptim, this.debug, this.funs, this.cfun, this.adjoint);
end        
