function rebuild_model_nested_events()
modelName = 'model_nested_events';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(1, 1, 5, 0, 0, 0, [], ['simulate_' modelName '.m'], 'model_nested_events', 'log10', 1, 1);
end
