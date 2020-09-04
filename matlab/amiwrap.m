function amiwrap( varargin )
    % AMIWRAP generates c++ mex files for the simulation of systems of differential equations via CVODES and IDAS.
    %
    % Parameters:
    %  varargin:
    %  modelname: specifies the name of the model which will be later used for the naming of the simulation file @type string
    %  symfun: specifies a function which executes model definition @type string.
    %  tdir: target directory where the simulation file should be placed @type string @default $AMICIDIR/models/modelname
    %  o2flag: boolean whether second order sensitivities should be enabled @type boolean @default false
    %
    % Return values:
    %  void

    matVer = ver('MATLAB');
    if(str2double(matVer.Version) >= 9.4)
        error('MATLAB R2018a or higher is currently not supported (see https://github.com/AMICI-dev/AMICI/issues/307)')
    end

    %%
    % check for MSVS
    if(~isempty(strfind(mex.getCompilerConfigurations('c++').Name,'Microsoft Windows')) || ~isempty(strfind(mex.getCompilerConfigurations('c++').Name,'Microsoft Visual')))
        warning('AMICI does not officially support Microsoft Visual Studio Compilers. If the compilation fails, we recommend using MinGW.')
    end

    %%
    % check inputs
    if(nargin<2)
        error('Must provide modelname and symfun.')
    end
    modelname = varargin{1}; % this is the target modelname
    if(~ischar(modelname))
        error('modelname must be a string.')
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
    % Display AMICI version
    disp(['amiwrap version ' getCommitHash(fileparts(mfilename('fullpath')))])

    %%
    % computations
    matlabRootPath=fileparts(mfilename('fullpath'));
    amiciRootPath=fileparts(matlabRootPath);

    addpath(genpath(fullfile(matlabRootPath,'auxiliary')));
    addpath(fullfile(matlabRootPath,'symbolic'));


    % try to load
    if(~isstruct(symfun))
        if(exist(symfun,'file')==2)
            model_hash = CalcMD5(which(symfun),'File');
        else
            model_hash = [];
        end
    else
        model_hash = [];
    end

    commit_hash = getCommitHash(amiciRootPath);

    if(~exist(fullfile(amiciRootPath,'models',modelname),'dir'))
        mkdir(fullfile(amiciRootPath,'models',modelname));
    end
    addpath(fullfile(amiciRootPath,'models',modelname));
    if(exist([commit_hash '_' model_hash '.mat'],'file')==2);
        load([commit_hash '_' model_hash '.mat']);
        % update modelname according to this function call
        model.updateModelName(modelname);
        % update wrap_path to this function call
        model.updateWrapPath(amiciRootPath);
    end

    if(~exist('model','var'))
        disp('Generating model struct ...')
        model = amimodel(symfun,modelname);


        if(~isempty(model_hash) && ~isempty(commit_hash))
            save(fullfile(amiciRootPath,'models',modelname,[commit_hash '_' model_hash]),'model')
        end
    end

    switch(o2flag)
        case 0
            o2string = [];
        case 1
            o2string = 'o2';
        case 2
            o2string = 'o2vec';
    end

    if(~isempty(o2string))
        o2_hash = CalcMD5(fullfile(matlabRootPath,'@amimodel',['augment' o2string '.m']),'File');
        try
            if(~exist(fullfile(amiciRootPath,'models',[modelname '_' o2string]),'dir'))
                mkdir(fullfile(amiciRootPath,'models',[modelname '_' o2string]));
            end
           addpath(fullfile(amiciRootPath,'models',[modelname '_' o2string]));
        end
        if(exist([commit_hash '_' model_hash '_' o2_hash '.mat'],'file')==2);
            load([commit_hash '_' model_hash '_' o2_hash '.mat']);
            % update modelname according to this function call
            modelo2.updateModelName([modelname '_' o2string]);
            % update wrap_path to this function call
            modelo2.updateWrapPath(amiciRootPath);
        end
        if(~exist('modelo2','var'))
            disp('Augmenting to second order ...')
            modelo2 = feval(['augment' o2string],model);


            if(~isempty(model_hash) && ~isempty(commit_hash))
                save(fullfile(amiciRootPath,'models',[modelname '_' o2string],[commit_hash '_' model_hash '_' o2_hash]),'modelo2')
            end
        end
    end

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
        clear(['ami_' modelname o2string]);
        movefile(fullfile(amiciRootPath,'models',modelname,['simulate_' modelname '.m']),fullfile(tdir,['simulate_' modelname '.m']));
        movefile(fullfile(amiciRootPath,'models',modelname,['ami_' modelname '.' mexext]),fullfile(tdir,['ami_' modelname '.' mexext]));
        % make files available in the path
        tmp = which(fullfile(tdir,['simulate_' modelname '.m']));
        tmp = which(fullfile(tdir,['ami_' modelname '.' mexext]));
        for fun = model.mfuns
            copyfile(fullfile(amiciRootPath,'models',modelname,[fun{1} '_' modelname '.m']),fullfile(tdir,[fun{1} '_' modelname '.m']));
            tmp = which(fullfile(tdir,[fun{1} '_' modelname '.m']));
        end
        % clear .m and .mex files from memory
        if(~isempty(o2string))
            movefile(fullfile(amiciRootPath,'models',[modelname '_' o2string],[ 'ami_' modelname '_' o2string '.' mexext]),fullfile(tdir,['ami_' modelname '_' o2string '.' mexext]));
            tmp = which(fullfile(tdir,['ami_' modelname '_' o2string '.' mexext]));
        end
    else
        addpath(fullfile(amiciRootPath,'models',modelname));
    end
    warning(warningreset);
end

