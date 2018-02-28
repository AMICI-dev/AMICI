function testModels()
    cd(fileparts(mfilename('fullpath')))
    addpath(genpath('cpputest'));
    wrapTestModels()
    cd(fileparts(mfilename('fullpath')))
    
    hdf5file = fullfile(fileparts(mfilename('fullpath')),'cpputest','expectedResults.h5');

    info = h5info(hdf5file);
    for imodel = 1:length(info.Groups)
        for itest = 1:length(info.Groups(imodel).Groups)
            [results,options,data] = readDataFromHDF5(info.Groups(imodel).Groups(itest),hdf5file);
            sol = getResults(info.Groups(imodel).Name(2:end),options,data);
            compareResults(sol,results);
        end 
    end
    
    
end

function [results,options,data] = readDataFromHDF5(group,hdf5file);
    data = [];
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

function sol = getResults(modelname,options,data)
    theta = options.theta;
    options = rmfield(options,'theta');
    kappa = options.kappa;
    options = rmfield(options,'kappa');
    t = options.ts;
    options = rmfield(options,'ts');
    options = rmfield(options,'nt');
    ami_options = amioption(options);
    if(~isempty(data))
        ami_data = amidata(data);
    else
        ami_data = [];
    end
    try
        sol = feval(['simulate_' modelname],t,theta,kappa,ami_data,ami_options);
    catch
        sol.status = -1;
    end
end

function compareResults(sol,results)
    if(results.status<0)
        assert(sol.status<0)
        return
    end

    for ifield = transpose(fieldnames(sol))
        if(strcmp(ifield{1},'diagnosis'))
            for jfield = transpose(fieldnames(sol.diagnosis))
                if(~ismember(jfield,{'J','newton_numsteps'})) % stuff that was previously missing/wrong
                    checkAgreement(sol.diagnosis.(jfield{1}),results.(jfield{1}),0,1);
                end
            end
        else
            if(~ismember(ifield,{'chi2','x0','sx0','ssigmay','s2x','s2y','s2sigmay'})) % stuff that was previously missing/wrong
                if(ismember(ifield,{'s2llh'}))
                    checkAgreement(sol.(ifield{1}),results.(ifield{1}),1e-4,1e-5)
                else
                    checkAgreement(sol.(ifield{1}),results.(ifield{1}))
                end
            end
        end
    end
end

function checkAgreement(actual,expected,atol,rtol)
   if(nargin<3)
       atol = 1e-10;
   end
   if(nargin<4)
       rtol = 1e-5;
   end
   if(~isempty(expected))
       assert(~isempty(actual));
       actual = actual(:);
       expected = expected(:);
       assert(all(isnan(actual)==isnan(expected)));
       actual = actual(~isnan(actual));
       expected = expected(~isnan(actual));
       assert(all(isinf(actual)==isinf(expected)));
       actual = actual(~isinf(actual));
       expected = expected(~isinf(actual));
       assert(all(abs(expected - actual) <= atol) || all(abs((expected - actual) ./ (rtol + expected)) <= rtol));
   end
end