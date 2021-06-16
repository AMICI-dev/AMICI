function rebuild_model_events()
modelName = 'model_events';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(3, 1, 4, 4, 2, 0, [], ['simulate_' modelName '.m'], 'model_events', 'log10', 1, 1);
end
