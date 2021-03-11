function rebuild_model_dirac()
modelName = 'model_dirac';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(2, 1, 4, 0, 0, 0, [], ['simulate_' modelName '.m'], 'model_dirac', 'log10', 1, 1);
end
