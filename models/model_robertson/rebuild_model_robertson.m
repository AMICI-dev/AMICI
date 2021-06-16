function rebuild_model_robertson()
modelName = 'model_robertson';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(3, 3, 3, 1, 0, 0, [], ['simulate_' modelName '.m'], 'model_robertson', 'log10', 1, 1);
end
