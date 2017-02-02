function amiwrap( varargin )
    % AMIWRAP generates c mex files for the simulation of systems of differential equations via CVODES and IDAS.
    %
    % Parameters:
    %  varargin:
    %  modelname: specifies the name of the model which will be later used for the naming of the simualation file @type string
    %  symfun: specifies a function which executes model defition see @ref definition for details @type string.
    %  tdir: target directory where the simulation file should be placed @type string @default $AMICIDIR/models/modelname
    %  o2flag: boolean whether second order sensitivities should be enabled @type boolean @default false
    %
    % Return values:
    %  void
    
    %%
    % check inputs
    if(nargin<2)
        error('Must provide modelname and symfun.')
    end
    modelname = varargin{1}; % this is the target modelname
    if(~ischar(modelname))
        error(' modelname must be a string')
    end
    symfun = varargin{2}; % this is the function which generates the symbolic struct
    if nargin > 2
        tdir = varargin{3};
    else
        tdir = pwd;
    end
    
    if nargin > 3
        o2flag = varargin{4};
        if(~ismember(o2flag,[0,1,2]))
            error('Parameter o2flag must have value 0, 1 or 2.');
        end
    else
        o2flag = false;
    end
    
    
    if(isempty(mex.getCompilerConfigurations('C')))
        error('No C compiler setup. Please install and configure with MATLAB')
    end
    if(~isempty(tdir))
        if(exist(tdir,'file') ~= 7)
            error('provided tdir is not a valid path')
        end
    end
    
    warningreset = warning;
    warning('off','symbolic:mupadmex:MuPADTextWarning')
    warning('off','MATLAB:dispatcher:nameConflict')
    warning('off','symbolic:sym:sym:DeprecateExpressions')
    warning('off','symbolic:generate:FunctionNotVerifiedToBeValid')
    %% 
    % computations
    
    % generate modelstruct
    disp('Generating model struct ...')
    model = amimodel(symfun,modelname);
    switch(o2flag)
        case 0 
            % do nothing
        case 1
            disp('Augmenting to second order ...')
            modelo2 = augmento2(model);
        case 2
            disp('Augmenting to second order ...')
            modelo2 = augmento2vec(model);
    end
    
    % do symbolic computations of modelstruct
    disp('Parsing model struct ...')
    model.parseModel();
    if(o2flag)
        modelo2.parseModel();
    end
    
    % generate C code out of symbolic computations
    disp('Generating C code ...')
    model.generateC();
    if(o2flag)
        modelo2.generateC();
    end
    
    % compile the previously generated C code
    disp('Compiling mex file ...')
    model.compileC();
    if(o2flag)
        modelo2.compileC();
    end
    
    % generate the matlab wrapper
    disp('Generating M code ...')
    if(o2flag)
        model.generateM(modelo2);
    else
        model.generateM([]);
    end
    
    if(~isempty(tdir))
        clear(['simulate_' modelname ]);
        clear(['ami_' modelname ]);
        clear(['ami_' modelname '_o2']);
        [odewrap_path,~,~]=fileparts(which('amiwrap.m'));
        movefile(fullfile(odewrap_path,'models',modelname,['simulate_' modelname '.m']),fullfile(tdir,['simulate_' modelname '.m']))
        movefile(fullfile(odewrap_path,'models',modelname,['ami_' modelname '.' mexext]),fullfile(tdir,['ami_' modelname '.' mexext]))
        for fun = model.mfuns
            copyfile(fullfile(odewrap_path,'models',modelname,[fun{1} '_' modelname '.m']),fullfile(tdir,[fun{1} '_' modelname '.m']))
        end
        % clear .m and .mex files from memory
        switch(o2flag)
            case 1
                movefile(fullfile(odewrap_path,'models',[modelname '_o2'],[ 'ami_' modelname '_o2.' mexext]),fullfile(tdir,['ami_' modelname '_o2.' mexext]))
            case 2
                movefile(fullfile(odewrap_path,'models',[modelname '_o2vec'],[ 'ami_' modelname '_o2vec.' mexext]),fullfile(tdir,['ami_' modelname '_o2vec.' mexext]))
        end

    end
    warning(warningreset);
end

