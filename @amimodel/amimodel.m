%
% @file amimodel
% @brief definition of amimodel class
%
classdef amimodel
    % amimodel is the object in which all model definitions are stored 
    
    properties ( GetAccess = 'public', SetAccess = 'private' )
        % symbolic definition struct @type struct
        sym@struct;
        % struct which stores information for which functions c code needs to be generated @type struct
        fun@struct;
        % struct which stores information for which functions c code needs
        % to be generated @type *amievent
        event@amievent;
        % name of the model @type string
        modelname@char;
        % struct that contains hash values for the symbolic model definitions @type struct
        HTable@struct;
        % default absolute tolerance @type double
        atol = 1e-8;
        % default relative tolerance @type double
        rtol = 1e-8;
        % default maximal number of integration steps @type int
        maxsteps = 1e4;
        % flag indicating whether debugging symbols should be compiled @type bool
        debug = false;
        % flag indicating whether adjoint sensitivities should be enabled @type bool
        adjoint = true;
        % flag indicating whether forward sensitivities should be enabled @type bool
        forward = true;
        % default initial time @type double
        t0 = 0;
        % type of wrapper (cvodes/idas) @type string
        wtype@char;
        % number of states @type int
        nx@double;
        % number of original states for second order sensitivities @type int
        nxtrue = 0;
        % number of observables @type int
        ny@double;
        % number of original observables for second order sensitivities @type int
        nytrue = 0;
        % number of parameters @type int
        np@double;
        % number of constants @type int
        nk@double;
        % number of events @type int
        nevent@double;
        % number of event outputs @type int
        nz@double;
        % flag for DAEs @type *int
        id@double;
        % upper Jacobian bandwidth @type int
        ubw@double;
        % lower Jacobian bandwidth @type int
        lbw@double;
        % number of nonzero entries in Jacobian @type int
        nnz@double;
        % dataindexes of sparse Jacobian @type *int
        sparseidx@double;
        % rowindexes of sparse Jacobian @type *int
        rowvals@double;
        % columnindexes of sparse Jacobian @type *int
        colptrs@double;
        % dataindexes of sparse Jacobian @type *int
        sparseidxB@double;
        % rowindexes of sparse Jacobian @type *int
        rowvalsB@double;
        % columnindexes of sparse Jacobian @type *int
        colptrsB@double;
        % cell array of functions to be compiled @type *cell
        funs@cell;
        % optimisation flag for compilation @type string
        coptim = '-O3';
        % default parametrisation @type string
        param = 'lin';
        % path to wrapper
        wrap_path@char;
        % flag to enforce recompilation of the model
        recompile = false;
        % storage for flags determining recompilation of individual
        % functions
        cfun@struct;
        % counter that allows enforcing of recompilation of models after
        % code changes
        compver = 2;
    end
    
    properties ( GetAccess = 'public', SetAccess = 'public' )
        % vector that maps outputs to events
        z2event@double;
    end
    
    methods
        function AM = amimodel(symfun,modelname)
            % constructor of the amimodel class. this function initializes the model object based on the provided
            % symfun and modelname
            %
            % Parameters:
            %  symfun: this is the string to the function which generates
            %  the modelstruct. You can also directly pass the struct here @type string
            %  modelname: name of the model @type string
            % 
            % Return values:
            %  AM: model definition object
            if(isa(symfun,'char'))
                model = eval(symfun);
            elseif(isa(symfun,'struct'))
                model = symfun;
            else
                error('invalid input symfun')
            end
            
            if(isfield(model,'sym'))
                AM.sym = model.sym;
            else
                error('symbolic definitions missing in struct returned by symfun')
            end
            
            props = properties(AM);
            
            for j = 1:length(props)
                if(~strcmp(props{j},'sym')) % we already checked for the sym field
                    if(isfield(model,props{j}))
                       AM.(props{j}) = model.(props{j});
                    end
                else
                    AM = AM.makeSyms();
                end
            end

            AM.modelname = modelname;
            % set path and create folder
            AM.wrap_path=fileparts(which('amiwrap.m'));
            if(~exist(fullfile(AM.wrap_path,'models'),'dir'))
                mkdir(fullfile(AM.wrap_path,'models'));
                mkdir(fullfile(AM.wrap_path,'models',AM.modelname));
            else
                if(~exist(fullfile(AM.wrap_path,'models',AM.modelname),'dir'))
                    mkdir(fullfile(AM.wrap_path,'models',AM.modelname))
                end
            end
            AM = AM.makeEvents();
            
            % check whether we have a DAE or ODE
            if(isfield(AM.sym,'M'))
                AM.wtype = 'iw'; % DAE
            else
                AM.wtype = 'cw'; % ODE
            end
        end
        
        this = parseModel(this)
        
        this = generateC(this)
        
        this = compileC(this)

        this = generateM(this,amimodelo2)
        
        this = getFun(this,HTable,funstr)
        
        this = makeEvents(this)
        
        this = makeSyms(this)
        
        [this,cflag] = checkDeps(this,HTable,deps)
        
        [this,HTable] = loadOldHashes(this) 
        
        [this] = augmento2(this)

    end
end

