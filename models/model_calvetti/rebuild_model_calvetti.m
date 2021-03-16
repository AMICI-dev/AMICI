function rebuild_model_calvetti()
modelName = 'model_calvetti';
amimodel.compileAndLinkModel(modelName, '', [], [], [], []);
amimodel.generateMatlabWrapper(6, 6, 0, 6, 0, 0, [], ['simulate_' modelName '.m'], 'model_calvetti', 'lin', 1, 1);
end
