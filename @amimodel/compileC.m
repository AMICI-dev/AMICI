function compileC(this)
    % compileC compiles the mex simulation file
    %
    % Return values:
    %  this: model definition object @type amimodel
    
    objectFileSuffix = '.o';
    if(ispc)
        objectFileSuffix = '.obj';
    end
    
    amiciSourcePath = fullfile(this.wrap_path,'src');
    modelSourceFolder = fullfile(this.wrap_path,'models',this.modelname);
    
    % compile flags
    COPT = ['COPTIMFLAGS=''' this.coptim ' -DNDEBUG'' CXXFLAGS=''$CXXFLAGS -std=c++0x'''];
    if(this.debug)
        DEBUG = ' -g CXXFLAGS=''$CXXFLAGS -Wall  -std=c++0x -Wno-unused-function -Wno-unused-variable'' ';
        COPT = ''; % no optimization with debug flags!
    else
        DEBUG = '';
    end
    
    %% Third party libraries
    [objectsstr, includesstr] = compileAMICIDependencies(this.wrap_path, objectFileSuffix, COPT, DEBUG);
    includesstr = strcat(includesstr,' -I"', modelSourceFolder, '"');
   
    %% Recompile AMICI base files if necessary
    % generate hash for file and append debug string if we have an md5
    % file, check this hash against the contained hash
    cppsrc = {'amici', 'symbolic_functions','spline', ...
        'edata','rdata','udata','tdata', ...
        'amici_interface_matlab', 'amici_misc', 'amici_model_functions'};
    objectArray = cellfun(@(x) [' "', fullfile(amiciSourcePath, x), objectFileSuffix, '"'], cppsrc, 'UniformOutput', false);
    objectsstr = [objectsstr, strjoin(objectArray, ' ')];
    sourcesForRecompile = cppsrc(cellfun(@(x) this.recompile || checkHash(fullfile(amiciSourcePath, x), objectFileSuffix, DEBUG), cppsrc));
    if(numel(sourcesForRecompile))
        fprintf('AMICI base files | ');
        sourceStr = '';
        for j = 1:numel(sourcesForRecompile)
            baseFilename = fullfile(amiciSourcePath, sourcesForRecompile{j});
            sourceStr  = [sourceStr, ' "', baseFilename, '.cpp"'];
        end
        eval(['mex ' DEBUG COPT ' -c -outdir ' amiciSourcePath ...
            includesstr ' ' sourceStr]);
        cellfun(@(x) updateFileHash(fullfile(amiciSourcePath, x), DEBUG), sourcesForRecompile);
    end
    
    %% Model-specific files
    for j=1:length(this.funs)
        baseFilename = fullfile(modelSourceFolder,[this.modelname '_' this.funs{j}]);
        recompile = this.recompile || checkHash(baseFilename,objectFileSuffix,DEBUG);
        this.cfun(1).(this.funs{j}) = recompile;
    end
    
    % flag dependencies for recompilation
    if(this.cfun(1).J)
        this.cfun(1).JBand = 1;
    end
    if(this.adjoint)
        if(this.cfun(1).JB)
            this.cfun(1).JBandB = 1;
        end
    end
    if(this.cfun(1).JSparse)
        this.cfun(1).sxdot = 1;
    end
    if(isfield(this.cfun(1),'dxdotdp'))
        if(this.cfun(1).dxdotdp)
            this.cfun(1).sxdot = 1;
            this.cfun(1).qBdot = 1;
        end
    end
    
    funsForRecompile = this.funs(structfun(@(x) logical(x), this.cfun(1)));
    if(numel(funsForRecompile))
        fprintf('ffuns | ');

        sources = cellfun(@(x) fullfile(modelSourceFolder,[this.modelname '_' x '.cpp']),funsForRecompile,'UniformOutput',false);
        sources = strjoin(sources,' ');
        
        eval(['mex ' DEBUG COPT ...
            ' -c -outdir ' modelSourceFolder ' ' ...
            sources ' ' ...
            includesstr ]);
        
        cellfun(@(x) updateFileHash(fullfile(modelSourceFolder,[this.modelname '_' x]), DEBUG),funsForRecompile,'UniformOutput',false);                
    end
    
    % append model object files
    for j=1:length(this.funs)
        objectsstr = strcat(objectsstr,...
            ' "',fullfile(modelSourceFolder, [this.modelname '_' this.funs{j} objectFileSuffix]),'"');
    end    
    
    % compile the wrapfunctions object
    fprintf('wrapfunctions | '); 
    eval(['mex ' DEBUG COPT ...
        ' -c -outdir ' modelSourceFolder ' ' ...
        fullfile(modelSourceFolder,'wrapfunctions.cpp') ' ' ...
        includesstr]);
    objectsstr = [objectsstr, ' "' fullfile(modelSourceFolder,['wrapfunctions' objectFileSuffix]) '"'];

    % now we have compiled everything model-specific, so we can replace hashes.mat to prevent recompilation
    try
        movefile(fullfile(modelSourceFolder,'hashes_new.mat'),...
        fullfile(modelSourceFolder,'hashes.mat'),'f');
    end
    
    %% Linking
    fprintf('linking | ');

    if(isunix)
        if(~ismac)
            CLIBS = 'CLIBS="-lrt -lmwblas"';
        else
            CLIBS = 'CLIBS="-lmwblas"';
        end
    else
        if(strcmp(mex.getCompilerConfigurations('c++').Name,'MinGW64 Compiler (C++)'))
            CLIBS = 'LD="g++"';
        else
            CLIBS = [];
        end
    end
    
    mexFilename = fullfile(modelSourceFolder,['ami_' this.modelname]);
    eval(['mex ' DEBUG ' ' COPT ' ' CLIBS ...
        ' -output ' mexFilename ' ' objectsstr])
end        
    
function hash = getFileHash(file)
    % getFileHash computed the md5hash of a given file
    %
    % Parameters:
    %  file: path of the file @type string
    %
    % Return values:
    %  hash: md5 hash of the provided file @type string
    hash = CalcMD5(file,'File','hex');
end    

function updateFileHash(baseFilename, DEBUG)
    hash = getFileHash([baseFilename '.cpp']);
    hash = [hash DEBUG];
    fid = fopen([baseFilename '_' mexext '.md5'],'w');
    fprintf(fid,hash);
    fclose(fid);
end
    
function recompile = checkHash(filestr,o_suffix,DEBUG)
    % checkHash checks whether filestr.cpp  has already been compiled as
    % filestr.o and whether the md5 hash of filestr.cpp matches the one in
    % filestr.md5
    %
    % Parameters:
    %  filestr: path of the file @type string
    %  o_suffix: OS specific suffix for compiled objects
    %  DEBUG: debug flag
    %
    % Return values:
    %  recompile: flag indicating whether we need to recompile filestr.cpp
    
    if(~exist([filestr o_suffix],'file'))
        % file does not exist, we need to recompile
        recompile = 1;
    else
        if(~exist([filestr '_' mexext '.md5'], 'file'))
            % hash does not exist, we need to recompile
            recompile = 1;
        else
            hash = getFileHash([filestr '.cpp']);
            hash = [hash DEBUG];
            fid = fopen([filestr '_' mexext '.md5']);
            tline = fgetl(fid);
            fclose(fid);
            if(~strcmp(tline,hash(1:end)))
                % file was updated, we need to recompile
                recompile = 1;
            else
                % everything is fine
                recompile = 0;
            end
        end
    end
end
