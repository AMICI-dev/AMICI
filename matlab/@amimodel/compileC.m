function compileC(this)
    % compileC compiles the mex simulation file
    %
    % Return values:
	%  void
	
    amimodel.compileAndLinkModel(this.modelname, fullfile(this.wrap_path,'models',this.modelname), this.coptim, this.debug, this.funs, this.cfun);
end        
