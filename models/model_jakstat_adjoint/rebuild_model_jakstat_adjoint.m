function rebuild_model_jakstat_adjoint()
modelName = 'model_jakstat_adjoint';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(9, 3, 17, 2, 0, 0, [], ['simulate_' modelName '.m'], 'model_jakstat_adjoint', 'log10', 1, 1);
end
