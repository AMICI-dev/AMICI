function rebuild_model_steadystate()
modelName = 'model_steadystate';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(3, 3, 5, 4, 0, 0, [], ['simulate_' modelName '.m'], 'model_steadystate', 'log10', 1, 1);
end
