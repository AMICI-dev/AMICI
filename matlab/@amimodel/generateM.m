function generateM(this, amimodelo2)
% generateM generates the matlab wrapper for the compiled C files.
%
% Parameters:
%  amimodelo2: this struct must contain all necessary symbolic
%  definitions for second order sensivities @type amimodel
%
% Return values:
%  void
wrapperFilename = fullfile(this.wrap_path,'models',this.modelname,['simulate_',this.modelname,'.m']);
this.generateMatlabWrapper(this.nx, this.ny, this.np, this.nk, this.nz, ...
    this.o2flag, amimodelo2, wrapperFilename, this.modelname, this.param, ...
    this.forward, this.adjoint)

%% Generation of the file which computes the Jacobian

for fun = this.mfuns
    if(isfield(this.fun,fun{1}))
        fprintf([fun{1} ' | ']);
        this.fun.(fun{1}).writeMcode(this);
    end
end


end

