function testModels()
    cd(fileparts(mfilename('fullpath')))
    addpath(genpath('cpputest'));
    wrapTestModels()
    cd(fileparts(mfilename('fullpath')))
    
    hdf5file = fullfile(fileparts(mfilename('fullpath')),'cpputest','expectedResults.h5');

    info = h5info(hdf5file);
    for imodel = 1:length(info.Groups)
        for itest = 1:length(info.Groups(imodel).Groups)
            [results,options,data,t,theta,kappa] = readDataFromHDF5(info.Groups(imodel).Groups(itest),hdf5file);
            sol = getResults(info.Groups(imodel).Name(2:end),options,data,t,theta,kappa);
            compareResults(sol,results);
        end 
    end
    
    
end

function [results,options,data,t,theta,kappa] = readDataFromHDF5(groups,hdf5file);
    data = [];
    t = [];
    theta = [];
    kappa = [];
    for igroup = 1:length(groups.Groups)
        group = groups.Groups(igroup);
        [~,name] = fileparts(group.Name);
        eval([name ' =  hdf2struct(group,''' group.Name ''',hdf5file);']);
    end
end

function s = hdf2struct(group,hdf5path,hdf5file)
    s = struct;
    if(~isempty(group.Attributes))
        for attr = {group.Attributes.Name};
            s.(attr{1}) = h5readatt(hdf5file,hdf5path,attr{1});
        end
    end
    if(~isempty(group.Datasets))
        for data = {group.Datasets.Name};
            s.(data{1}) = h5read(hdf5file,[hdf5path '/' data{1}]);
        end
    end
    for igroup = 1:length(group.Groups);
        [~,name] = fileparts(group.Groups(igroup).Name);
        s.(name) = hdf2struct(group.Groups(igroup),[hdf5path '/' name],hdf5file);
    end
end

function sol = getResults(modelname,options,data,t,theta,kappa)
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
                if(~ismember(jfield{1},{'newton_time'}))
                    checkAgreement(sol.diagnosis,results.diagnosis,jfield{1},0,1);
                end
            end
        else
            if(ismember(ifield,{'s2llh', 's2x', 's2y', 's2z', 'sx'}))
                checkAgreement(sol,results,ifield{1},1e-4,1e-3)
            else
                checkAgreement(sol,results,ifield{1})
            end
        end
    end
end

function checkAgreement(sol,results,fieldname,atol,rtol)
    if(~isfield(results,fieldname))
        assert(isempty(sol.(fieldname)))
        return
    end
    expected = results.(fieldname);
    actual = sol.(fieldname);
    if(nargin<4)
        atol = 1e-10;
    end
    if(nargin<5)
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
        assert(all(abs(expected - actual) <= atol) || all(abs((expected - actual) ./ (rtol + abs(expected))) <= rtol));
    end
end