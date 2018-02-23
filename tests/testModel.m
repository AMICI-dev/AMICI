function testModels()
%     wrapTestModels()

    hdf5file = './cpputest/expectedResults.h5'

    info = h5info(hdf5file)
    for imodel = 1:length(info.Groups)
        for itest = 1:length(info.Groups(imodel).Groups)
            [results,options] = readDataFromHDF5(info.Groups(imodel).Groups(itest),hdf5file);
            sol = getResults(info.Groups(imodel).Name(2:end),options)
            compareResults(sol,results);
        end 
    end
    
    
end

function [results,options] = readDataFromHDF5(group,hdf5file);
    for idataset = 1:length(group.Datasets)
        dataset = group.Datasets(idataset);
        eval([dataset.Name ' =  hdf2struct(dataset,''' group.Name '/' dataset.Name ''',hdf5file);']);
    end
end

function s = hdf2struct(dataset,hdf5path,hdf5file)
    s = struct;
    for attr = {dataset.Attributes.Name};
        s.(attr{1}) = h5readatt(hdf5file,hdf5path,attr{1});
    end
end

function sol = getResults(modelname,options)
    theta = options.theta;
    options = rmfield(options,'theta')
    kappa = options.kappa;
    options = rmfield(options,'kappa')
    t = options.ts;
    options = rmfield(options,'ts')
    options = rmfield(options,'nt')
    ami_options = amioption(options)
    sol = feval(['simulate_' modelname],t,theta,kappa,[],ami_options);
end

function compareResults(sol,results)

    for ifield = transpose(fieldnames(sol))
        if(strcmp(ifield{1},'diagnosis'))
            for jfield = transpose(fieldnames(sol.diagnosis))
                checkAgreement(sol.diagnosis.(jfield{1}),results(jfield{1}));
            end
        else
            checkAgreement(sol.(ifield{1}),results.(ifield{1}))
        end
    end
end

function checkAgreement(actual,expected)
   atol = 1e-5;
   rtol = 1e-5;
   if(~isempty(expected))
       assert(~isempty(actual))
       actual = actual(:);
       expected = expected(:);
       assert(isnan(actual)==isnan(expected))
       actual = actual(~isnan(actual));
       expected = expected(~isnan(actual));
       assert(all(abs(expected - actual) <= atol) || all(abs((expected - actual) / (rtol + expected)) <= rtol))
   end
end