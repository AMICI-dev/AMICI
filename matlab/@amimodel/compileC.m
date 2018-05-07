function compileC(this)
    % compileC compiles the mex simulation file
    %
    % Return values:
	%  void
	
    amimodel.compileAndLinkModel(this.modelname, this.wrap_path, this.recompile, this.coptim, this.debug, this.funs, this.cfun, this.adjoint);
end        
