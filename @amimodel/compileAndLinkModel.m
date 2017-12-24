function compileAndLinkModel(modelname, wrap_path, recompile, coptim, debug, funs, cfun, adjoint)
    % compileAndLinkModel compiles the mex simulation file.
    % It does not check if the model files have changed since generating 
    % C++ code or whether all files are still present. 
    % Use only if you know what you are doing. The safer alternative is 
    % rerunning amiwrap().
    % 
    % Parameters:
    %  modelname: name of the model as specified for amiwrap()
    %  wrap_path: AMICI path
    %  recompile: flag indicating whether all source files should be
    %  recompiled
    %  coptim: optimization flags
    %  debug: enable debugging
    %  funs: array with names of the model functions, will be guessed 
    %   from source files if left empty 
    %  cfun: struct indicating which files should be recompiled
    %  adjoint: flag indicating whether adjoint sensitivies are enabled
    %
    % Return values:
    %  void
    
    % if no list provided, try to determine relevant files from model
    % folder
    if(isempty(funs))
        ls = dir(fullfile(wrap_path, 'models', modelname, [modelname '_*.cpp']));
        ls = {ls.name};
        % extract funs from filename (strip of modelname_ and .cpp
        funs = cellfun(@(x) x((length(modelname)+2):(length(x)-4)), ls, 'UniformOutput', false);
    end
    
    objectFileSuffix = '.o';
    if(ispc)
        objectFileSuffix = '.obj';
    end
    
    amiciSourcePath = fullfile(wrap_path,'src');
    modelSourceFolder = fullfile(wrap_path,'models',modelname);
    
    % compile flags
    COPT = ['COPTIMFLAGS=''' coptim ' -DNDEBUG'' CXXFLAGS=''$CXXFLAGS -std=c++0x'''];
    if(debug)
        DEBUG = ' -g CXXFLAGS=''$CXXFLAGS -Wall  -std=c++0x -Wno-unused-function -Wno-unused-variable'' ';
        COPT = ''; % no optimization with debug flags!
    else
        DEBUG = '';
    end
    
    %% Third party libraries
    [objectsstr, includesstr] = compileAMICIDependencies(wrap_path, objectFileSuffix, COPT, DEBUG);
    includesstr = strcat(includesstr,' -I"', modelSourceFolder, '"');
   
    %% Recompile AMICI base files if necessary
    [objectStrAmici, recompile ] = compileAmiciBase(amiciSourcePath, objectFileSuffix, includesstr, recompile, DEBUG, COPT);
    objectsstr = [objectsstr, objectStrAmici];

    %% Model-specific files
    for j=1:length(funs)
        baseFilename = fullfile(modelSourceFolder,[modelname '_' funs{j}]);
        recompile = recompile || sourceNeedsRecompilation(baseFilename,objectFileSuffix,DEBUG);
        cfun(1).(funs{j}) = recompile;
    end
    
    funsForRecompile = {};
    
    % flag dependencies for recompilation
    if(~isempty(cfun))
        if(isfield('J',cfun(1)))
            if(cfun(1).J)
                if(ismember('JBand',funs))
                    cfun(1).JBand = 1;
                end
            end
        end
        
        if(isfield('JB',cfun(1)))
            if(cfun(1).JB)
                if(ismember('JBandB',funs))
                    cfun(1).JBandB = 1;
                end
            end
        end
        
        if(isfield('JSparse',cfun(1)))
            if(cfun(1).JSparse)
                if(ismember('sxdot',funs))
                    cfun(1).sxdot = 1;
                end
            end
        end
        if(isfield('dxdotdp',cfun(1)))
            if(isfield(cfun(1),'dxdotdp'))
                if(cfun(1).dxdotdp)
                    if(ismember('sxdot',funs))
                        cfun(1).sxdot = 1;
                    end
                    if(ismember('qBdot',funs))
                        cfun(1).qBdot = 1;
                    end
                end
            end
        end
        funsForRecompile = funs(structfun(@(x) logical(x), cfun(1)));
    end
    
    if(numel(funsForRecompile))
        fprintf('ffuns | ');

        sources = cellfun(@(x) fullfile(modelSourceFolder,[modelname '_' x '.cpp']),funsForRecompile,'UniformOutput',false);
        sources = strjoin(sources,' ');
        
        eval(['mex ' DEBUG COPT ...
            ' -c -outdir ' modelSourceFolder ' ' ...
            sources ' ' ...
            includesstr ]);
        
        cellfun(@(x) updateFileHashSource(fullfile(modelSourceFolder,[modelname '_' x]), DEBUG),funsForRecompile,'UniformOutput',false);                
    end
    
    % append model object files
    for j=1:length(funs)
        objectsstr = strcat(objectsstr,...
            ' "',fullfile(modelSourceFolder, [modelname '_' funs{j} objectFileSuffix]),'"');
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
            CLIBS = '-lmwblas';
        end
    end
    
    mexFilename = fullfile(modelSourceFolder,['ami_' modelname]);
    eval(['mex ' DEBUG ' ' COPT ' ' CLIBS ...
        ' -output ' mexFilename ' ' objectsstr])
end        
    
function [objectStrAmici, recompile ] = compileAmiciBase(amiciSourcePath, objectFileSuffix, includesstr, recompile, DEBUG, COPT)
    % generate hash for file and append debug string if we have an md5
    % file, check this hash against the contained hash
    cppsrc = {'amici', 'symbolic_functions','spline', ...
        'edata','rdata','udata','tdata', ...
        'amici_interface_matlab', 'amici_misc', ...
        'amici_solver', 'amici_solver_cvodes', 'amici_solver_idas', ...
        'amici_model', 'returndata_matlab', ...
        'forwardproblem', 'steadystateproblem', 'backwardproblem', 'newton_solver'};
    % to be safe, recompile everything if headers have changed. otherwise
    % would need to check the full include hierarchy
    recompile = recompile || headersHaveChanged([amiciSourcePath '/../include/'], DEBUG);
    objectArray = cellfun(@(x) [' "', fullfile(amiciSourcePath, x), objectFileSuffix, '"'], cppsrc, 'UniformOutput', false);
    objectStrAmici = strjoin(objectArray, ' ');
    sourcesForRecompile = cppsrc(cellfun(@(x) recompile || sourceNeedsRecompilation(fullfile(amiciSourcePath, x), objectFileSuffix, DEBUG), cppsrc));
    if(numel(sourcesForRecompile))
        fprintf('AMICI base files | ');
        sourceStr = '';
        for j = 1:numel(sourcesForRecompile)
            baseFilename = fullfile(amiciSourcePath, sourcesForRecompile{j});
            sourceStr  = [sourceStr, ' "', baseFilename, '.cpp"'];
        end
        eval(['mex ' DEBUG COPT ' -c -outdir ' amiciSourcePath ...
            includesstr ' ' sourceStr]);
        cellfun(@(x) updateFileHashSource(fullfile(amiciSourcePath, x), DEBUG), sourcesForRecompile);
        updateHeaderFileHashes([amiciSourcePath '/../include/'], DEBUG);
    end
    
end

function headersChanged = headersHaveChanged(includePath, DEBUG)
    list = dir([includePath '/*.h']);
    headersChanged = false;
    for file = {list.name}
        headersChanged = headerFileChanged(fullfile(includePath, file{:}), DEBUG);
        if(headersChanged)
            break;
        end            
    end
end

function updateHeaderFileHashes(includePath, DEBUG)
    list = dir([includePath '/*.h']);
    for file = {list.name}
        updateFileHash(fullfile(includePath, file{:}), DEBUG);
    end
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

function updateFileHashSource(baseFilename, DEBUG)
    if(exist([baseFilename '.cpp'], 'file'))
        updateFileHash([baseFilename '.cpp'], DEBUG);
    end
end


function updateFileHash(filename, DEBUG)
    hash = getFileHash(filename);
    hash = [hash DEBUG];
    fid = fopen([filename '_' mexext '.md5'],'w');
    fprintf(fid,hash);
    fclose(fid);
end

function headerChanged = headerFileChanged(filename, DEBUG)
    hashFileSufffix = ['_' mexext '.md5'];
    
    if(~exist(filename, 'file') && exist([filename, hashFileSufffix], 'file'))
        % header file has been removed, recompile and remove stray
        % hash file
        headerChanged = true;
        delete([filename, hashFileSufffix])      
    elseif(~exist([filename, hashFileSufffix], 'file'))
        % there exists a header file but no corresponding hash file
        headerChanged = true;
    else
        % hash file exist, did hash change?
        headerChanged = hashHasChanged(filename, [filename hashFileSufffix], DEBUG);
    end
end

function recompile = sourceNeedsRecompilation(filestr, o_suffix,DEBUG)
    % sourceNeedsRecompilation checks whether filestr.cpp  has already been 
    % compiled as filestr.o and whether the md5 hash of filestr.cpp matches 
    % the one in filestr.md5
    %
    % Parameters:
    %  * filestr: path of the file @type string
    %  * o_suffix: OS specific suffix for compiled objects
    %  * DEBUG: debug flag
    %
    % Return values:
    %  recompile: flag indicating whether we need to recompile filestr.cpp
    
    if(~exist([filestr o_suffix],'file'))
        % object file does not exist, we need to recompile
        recompile = 1;
    else
        sourceFilename = [filestr '.cpp'];
        hashFileName = [sourceFilename '_' mexext '.md5'];
        if(~exist(hashFileName, 'file'))
            % source file hash does not exist, we need to recompile
            recompile = 1;
        else
            % hash files exists, did they change?
            recompile = hashHasChanged(sourceFilename, hashFileName, DEBUG);
        end
    end
end

function hasChanged = hashHasChanged(sourceFilename, hashFilename, hashSuffix)
    % checkHash checks whether the given file matches the saved hash
    % 
    % Parameters:
    %  * sourceFilename: the file to hash (has to exist)
    %  * hashFilename: the file where the hash is saved (has to exist)
    %  * hashSuffix: string to be appended to the hash
    % Return values:
    %  * hasChanged: true if file and saved hash do not match

    hash = getFileHash(sourceFilename);
    hash = [hash hashSuffix];
    fid = fopen(hashFilename);
    tline = fgetl(fid);
    fclose(fid);
    hasChanged = ~strcmp(tline, hash(1:end));
end
