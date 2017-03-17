%
% @file amised
% @brief definition of amised class
%
classdef amised < handle
    % AMISED is a container for SED-ML objects
    
    properties ( GetAccess = 'public', SetAccess = 'private' )
        % amimodel from the specified model
        model = struct('event',[],'sym',[]);
        % cell array of model identifiers
        modelname = {};
        % stores the struct tree from the xml definition
        sedml = struct.empty();
        % count the number of outputs per model
        outputcount = [];
        % indexes for dataGenerators
        varidx = [];
        % symbolic expressions for variables
        varsym = sym([]);
        % symbolic expressions for data
        datasym = sym([]);
        
    end
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
 
    end
    
    methods
        function ASED = amised(sedname)
            %amised reads in an SEDML document using the JAVA binding of 
            % of libSEDML
            %
            % Parameters:
            %  sedname: name/path of the SEDML document
            %
            % Return values:
            %  ASED: amised object which contains all the information from
            %        the SEDML document
            
            % get models
            for imodel = 1:length(ASED.sedml.listOfModels.model)
                % get the model sbml
                % check if this is a biomodels model
                if(length(ASED.sedml.listOfModels.model{imodel}.Attributes.source>=23))
                    if(strcmp(ASED.sedml.listOfModels.model{imodel}.Attributes.source,ASED.modelname))
                        ASED.model(imodel) = ASED.model(find(strcmp(ASED.sedml.listOfModels.model{imodel}.Attributes.source,ASED.modelname)));
                        ASED.modelname{imodel} = ASED.sedml.listOfModels.model{imodel}.Attributes.id;
                        ASED.model(imodel).sym.y = sym([]);
                        ASED.outputcount(imodel) = 0;
                        continue
                    end
                    if(strcmp(ASED.sedml.listOfModels.model{imodel}.Attributes.source(1:23),'urn:miriam:biomodels.db'))
                        modelxml = websave([ASED.sedml.listOfModels.model{imodel}.Attributes.source(25:end) '.xml'],['http://www.ebi.ac.uk/biomodels-main/download?mid=' ASED.sedml.listOfModels.model{imodel}.Attributes.source(25:end)]);
                    else
                        modelxml = ASED.sedml.listOfModels.model{imodel}.Attributes.source;
                    end
                else
                    modelxml = ASED.sedml.listOfModels.model{imodel}.Attributes.source;
                end
                modelxml = strrep(modelxml,'.xml','');
                % convert model sbml to amici
                SBML2AMICI(modelxml,ASED.sedml.listOfModels.model{imodel}.Attributes.id);
                eval(['model = ' ASED.sedml.listOfModels.model{imodel}.Attributes.id '_syms();'])
                if(~isfield(model,'event'))
                    model.event = [];
                end
                ASED.model(imodel) = model;
                % clear output;
                ASED.model(imodel).sym.y = sym([]);
                ASED.outputcount(imodel) = 0;
                ASED.modelname{imodel} = ASED.sedml.listOfModels.model{imodel}.Attributes.id;
            end
            % apply changes
            if(isfield(ASED.sedml,'listOfChanges'))
                %TBD apply changes
            end
            % construct outputs
            for idata = 1:length(ASED.sedml.listOfDataGenerators.dataGenerator)
                if(length(ASED.sedml.listOfDataGenerators.dataGenerator)>1)
                    dataGenerator = ASED.sedml.listOfDataGenerators.dataGenerator{idata};
                else
                    dataGenerator = ASED.sedml.listOfDataGenerators.dataGenerator;
                end
                tasks = [ASED.sedml.listOfTasks.task{:}];
                tasksatts = [tasks.Attributes];
                taskids = {tasksatts.id};
                for ivar = 1:length(dataGenerator.listOfVariables.variable)
                    if(length(dataGenerator.listOfVariables.variable)>1)
                        variable = dataGenerator.listOfVariables.variable{ivar};
                    else
                        variable = dataGenerator.listOfVariables.variable;
                    end
                    task_id = find(strcmp(variable.Attributes.taskReference,taskids));
                    if(isempty(task_id))
                        error(['Invalid taskReference in dataGenerator ' num2str(idata) ': ' variable.Attributes.taskReference]);
                    end
                    model_idx = find(strcmp(ASED.sedml.listOfTasks.task{task_id}.Attributes.modelReference,ASED.modelname));
                    if(isempty(model_idx))
                        error(['Invalid modelReference in task ' num2str(task_id) ': ' ASED.sedml.listOfTasks.task{task_id}.Attributes.modelReference]);
                    end
                    if(isfield(variable.Attributes,'symbol'))
                        if(strcmp(variable.Attributes.symbol,'urn:sedml:symbol:time'))
                            ASED.model(model_idx).sym.y(ASED.outputcount(model_idx)+1) = sym('t');
                        end
                    end
                    if(isfield(variable.Attributes,'target'))
                        if(strfind(variable.Attributes.target,'/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species'))
                            ASED.model(model_idx).sym.y(ASED.outputcount(model_idx)+1) = sym(variable.Attributes.target(60:end-2));
                        end
                    end
                    ASED.outputcount(model_idx) = ASED.outputcount(model_idx)+1;
                    ASED.varidx(idata,ivar,:) = [model_idx,ASED.outputcount(model_idx)];
                    ASED.varsym(idata,ivar) = sym(variable.Attributes.id);
                end
                ASED.datasym(idata) = sym(variable.Attributes.id);
                
            end

            
        end
    end
end

