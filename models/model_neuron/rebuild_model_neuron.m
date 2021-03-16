function rebuild_model_neuron()
modelName = 'model_neuron';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(2, 1, 4, 2, 1, 0, [], ['simulate_' modelName '.m'], 'model_neuron', 'log10', 1, 1);
end
