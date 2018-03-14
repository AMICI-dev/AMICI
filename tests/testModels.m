function testModels()
    
    % disable specific warnings for these tests, some tests are supposed
    % to produce warnings
    warningreset = warning;
    warning('off','AMICI:mex:simulation')
    warning('off','AMICI:mex:CVODES:CVode:TOO_MUCH_WORK')
    
    ignoredTests = {'/model_jakstat_adjoint/sensiadjointemptysensind', ...
                    '/model_jakstat_adjoint/sensiforwardemptysensind'};
    
    cd(fileparts(mfilename('fullpath')))
    addpath(genpath('cpputest'));
    wrapTestModels()
    cd(fileparts(mfilename('fullpath')))
    
    hdf5file = fullfile(fileparts(mfilename('fullpath')),'cpputest','expectedResults.h5');

    info = h5info(hdf5file);
    for imodel = 1:length(info.Groups)
        if(~isempty(regexp(info.Groups(imodel).Name(2:end),'^model_neuron')))
            model_atol = 1e-9;
            model_rtol = 1e-4;
        else
            model_atol = 1e-10;
            model_rtol = 1e-5;
        end
        for itest = 1:length(info.Groups(imodel).Groups)
            if(ismember(info.Groups(imodel).Groups(itest).Name, ignoredTests))
                continue
            end
            
            [results,options,data,t,theta,kappa] = readDataFromHDF5(info.Groups(imodel).Groups(itest),hdf5file);
            sol = getResults(info.Groups(imodel).Name(2:end),options,data,t,theta,kappa);
            compareResults(sol,results);
        end
    end
    
    warning(warningreset);
    
    %% begin nested functions
    
    function sol = getResults(modelname,options,data,t,theta,kappa)
        theta = options.theta;
        options = rmfield(options,'theta');
        kappa = options.kappa;
        options = rmfield(options,'kappa');
        t = options.ts;
        options = rmfield(options,'ts');
        ami_options = amioption(options);
        if(~isempty(data))
            ami_data = amidata(data);
        else
            ami_data = [];
        end
        sol = feval(['simulate_' modelname],t,theta,kappa,ami_data,ami_options);
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
                if(~ismember(ifield,{'chi2'}))
                    if(~isempty(sol.s2llh))
                        if(ismember(ifield,{'s2llh'}))
                            checkAgreement(sol,results,ifield{1},1e-4,1e-3)
                        elseif(ismember(ifield,{'x','sx','s2x','y','sy','s2y','z','sz','s2z','rz','srz','s2rz','sigmay','ssigmay','s2sigmay','sigmaz','ssigmaz','s2sigmaz','x0','sx0'}))
                            % TODO: reimplement this as soon as #264 is done
                        else
                            checkAgreement(sol,results,ifield{1})
                        end
                    else
                        checkAgreement(sol,results,ifield{1})
                    end
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
            atol = model_atol;
        end
        if(nargin<5)
            rtol = model_rtol;
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
        if(isfield(options, 'sens_ind'))
            % adapt to base 1 indexing
            options.sens_ind = options.sens_ind + 1;
        end
        if(~isempty(data))
            if(length(data.t)~=size(data.Y,1)) % no idea why this goes wrong only _sometimes_
                data.Y = transpose(data.Y);
                data.Sigma_Y = transpose(data.Sigma_Y);
            end
            if(size(data.Z,2)>0) % don't ask ...
                data.Z = transpose(data.Z);
                data.Sigma_Z = transpose(data.Sigma_Z);
            end
        end
    end
    
    
    
    function matlab = cpp2matlab(cpp)
        dims = size(cpp);
        if(sum(dims>1)>1 || length(dims)>2)
            switch(length(dims))
                case 3
                    matlab = permute(cpp,[3,1,2]);
                case 2
                    matlab = transpose(cpp);
                otherwise
                    matlab = cpp;
            end
        else
            matlab = cpp;
        end
        matlab = double(matlab);
    end
    
    function s = hdf2struct(group,hdf5path,hdf5file)
        s = struct;
        if(~isempty(group.Attributes))
            for attr = {group.Attributes.Name};
                s.(attr{1}) = double(h5readatt(hdf5file,hdf5path,attr{1}));
            end
        end
        if(~isempty(group.Datasets))
            for data = {group.Datasets.Name};
                s.(data{1}) = cpp2matlab(h5read(hdf5file,[hdf5path '/' data{1}]));
            end
        end
        for igroup = 1:length(group.Groups);
            [~,name] = fileparts(group.Groups(igroup).Name);
            s.(name) = hdf2struct(group.Groups(igroup),[hdf5path '/' name],hdf5file);
        end
    end
end
