function  SBML2AMICI( filename, modelname )
    % SBML2AMICI generates AMICI model definition files from SBML.
    %
    % Parameters:
    %  filename: name of the SBML file (withouth extension)
    %  modelname: name of the model, this will define the name of the
    %  output file (default: input filename)
    %
    % Return values:
    %  void
    if(nargin<2)
        modelname = filename;
    end
    wrap_path=fileparts(mfilename('fullpath'));
    if(~exist(fullfile(wrap_path,'SBMLimporter'),'dir'))
        error('SBMLimporter is not available, try running `git submodule update --init`')
    end
    addpath(fullfile(wrap_path,'SBMLimporter'));
    ODE = SBMLode(filename);
    ODE.writeAMICI(modelname);
    pnom = ODE.pnom;
    knom = ODE.knom;
    vnom = ODE.volume;
    kvnom = ODE.kvolume;
    save([modelname '_pnom.mat'],'pnom');
    save([modelname '_knom.mat'],'knom');
    save([modelname '_vnom.mat'],'vnom');
    save([modelname '_kvnom.mat'],'kvnom');
end

