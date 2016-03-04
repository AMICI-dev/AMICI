function  SBML2AMICI( filename, modelname )
    % SBML2AMICI generates AMICI model definition files from SBML.
    %
    % Parameters:
    %  filename: name of the SBML file (withouth extension)
    %  modelname: name of the model, this will define the name of the
    %  output file
    %
    % Return values:
    %  void
    
    ODE = SBMLode(filename);
    ODE.writeAMICI(modelname);
    
    
end

