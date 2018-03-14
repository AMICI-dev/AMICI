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
    
    if(debug)
        objectFolder = fullfile(wrap_path,'models',mexext,'debug');
    else
        objectFolder = fullfile(wrap_path,'models',mexext,'release');
    end
    
    %% Third party libraries
    [objectsstr, includesstr] = compileAMICIDependencies(wrap_path, objectFolder, objectFileSuffix, COPT, DEBUG);
    includesstr = strcat(includesstr,' -I"', modelSourceFolder, '"');
   
    %% Recompile AMICI base files if necessary
    [objectStrAmici] = compileAmiciBase(wrap_path, objectFolder, objectFileSuffix, includesstr, DEBUG, COPT);
    objectsstr = [objectsstr, objectStrAmici];

    %% Model-specific files
    for j=1:length(funs)
        baseFilename = fullfile(modelSourceFolder,[modelname '_' funs{j}]);
        cfun(1).(funs{j}) = sourceNeedsRecompilation(baseFilename,objectFileSuffix,DEBUG);
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
        filename = fullfile(modelSourceFolder, [modelname '_' funs{j} objectFileSuffix]);
        if(exist(filename,'file'))
            objectsstr = strcat(objectsstr,...
                ' "',filename,'"');
        end
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
            CLIBS = 'CLIBS="-lrt -lmwblas -ldl"';
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
    
function [objectStrAmici] = compileAmiciBase(wrap_path, objectFolder, objectFileSuffix, includesstr, DEBUG, COPT)
    % generate hash for file and append debug string if we have an md5
    % file, check this hash against the contained hash
    cppsrc = {'amici', 'symbolic_functions','spline', ...
        'edata','rdata', ...
        'interface_matlab', 'misc', ...
        'solver', 'solver_cvodes', 'solver_idas', ...
        'model', 'model_ode', 'model_dae', 'returndata_matlab', ...
        'forwardproblem', 'steadystateproblem', 'backwardproblem', 'newton_solver'};
    % to be safe, recompile everything if headers have changed. otherwise
    % would need to check the full include hierarchy
    amiciIncludePath = fullfile(wrap_path,'include','amici');
    amiciSourcePath = fullfile(wrap_path,'src');
    recompile = headersHaveChanged(amiciIncludePath,objectFolder);
    objectArray = cellfun(@(x) [' "', fullfile(objectFolder, x), objectFileSuffix, '"'], cppsrc, 'UniformOutput', false);
    objectStrAmici = strjoin(objectArray, ' ');
    sourcesForRecompile = cppsrc(cellfun(@(x) recompile || sourceNeedsRecompilation(amiciSourcePath, objectFolder, x, objectFileSuffix), cppsrc));
    if(numel(sourcesForRecompile))
        fprintf('AMICI base files | ');
        sourceStr = '';
        for j = 1:numel(sourcesForRecompile)
            baseFilename = fullfile(amiciSourcePath, sourcesForRecompile{j});
            sourceStr  = [sourceStr, ' "', baseFilename, '.cpp"'];
        end
        eval(['mex ' DEBUG COPT ' -c -outdir ' objectFolder ...
            includesstr ' ' sourceStr]);
        cellfun(@(x) updateFileHashSource(amiciSourcePath, objectFolder, x), sourcesForRecompile);
        updateHeaderFileHashes(amiciIncludePath, objectFolder);
    end
    
end

function headersChanged = headersHaveChanged(includePath, objectFolder)
    list = dir([includePath '/*.h']);
    headersChanged = false;
    for file = {list.name}
        headersChanged = headerFileChanged(includePath, objectFolder, file{:});
        if(headersChanged)
            break;
        end            
    end
end

function updateHeaderFileHashes(includePath, objectFolder)
    list = dir([includePath '/*.h']);
    for file = {list.name}
        updateFileHash(includePath, objectFolder, file{:});
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

function updateFileHashSource(sourceFolder,objectFolder,baseFilename)
    fileName = [baseFilename '.cpp'];
    if(exist(fullfile(sourceFolder,fileName), 'file'))
        updateFileHash(sourceFolder,objectFolder,fileName);
    end
end


function updateFileHash(fileFolder,hashFolder,filename)
    hash = getFileHash(fullfile(fileFolder,filename));
    fid = fopen(fullfile(hashFolder,[filename '.md5']),'w');
    fprintf(fid,hash);
    fclose(fid);
end

function headerChanged = headerFileChanged(includePath, objectFolder, fileName)
    % headerFileChanged checks whether fileName.h  has changed since last compilation
    %
    % Parameters:
    %  * includePath: path to directory containing header file @type string
    %  * objectFolder: path to directory containing compiled md5 file @type string
    %  * fileName: name of header @type string
    %
    % Return values:
    %  headerChanged: flag indicating whether we need to recompile filestr.cpp
    hashFileSufffix = ['.md5'];
    headerFile = fullfile(includePath,fileName);
    hashFile = fullfile(objectFolder,[fileName hashFileSufffix]);
    if(~exist(headerFile , 'file') && exist(hashFile, 'file'))
        % header file has been removed, recompile and remove stray
        % hash file
        headerChanged = true;
        delete(hashFile);
    elseif(~exist(hashFile, 'file'))
        % there exists a header file but no corresponding hash file
        headerChanged = true;
    else
        % hash file exist, did hash change?
        headerChanged = hashHasChanged(headerFile, hashFile);
    end
end

function recompile = sourceNeedsRecompilation(amiciSourcePath, objectFolder, fileName, o_suffix)
    % sourceNeedsRecompilation checks whether fileName.cpp  has already been 
    % compiled as fileName.o and whether the md5 hash of fileName.cpp matches 
    % the one in fileName.md5
    %
    % Parameters:
    %  * amiciSourcePath: path to directory containing source file @type string
    %  * objectFolder: path to directory containing compiled object and md5 file @type string
    %  * fileName: name of source @type string
    %  * o_suffix: OS specific suffix for compiled objects @type string
    %
    % Return values:
    %  recompile: flag indicating whether we need to recompile filestr.cpp
    
    sourceFile = fullfile(amiciSourcePath,[fileName '.cpp']);
    
    if(~exist(sourceFile,'file'))
        % cpp does not exist, we don't need to compile :)
        recompile = 0;
    elseif(~exist(fullfile(objectFolder,[fileName o_suffix]),'file'))
        % object file does not exist, we need to recompile
        recompile = 1;
    else
        hashFile = fullfile(objectFolder,[fileName '.cpp.md5']);
        if(~exist(hashFile, 'file'))
            % source file hash does not exist, we need to recompile
            recompile = 1;
        else
            % hash files exists, did they change?
            recompile = hashHasChanged(sourceFile, hashFile);
        end
    end
end

function hasChanged = hashHasChanged(sourceFilename, hashFilename)
    % checkHash checks whether the given file matches the saved hash
    % 
    % Parameters:
    %  * sourceFilename: the file to hash (has to exist)
    %  * hashFilename: the file where the hash is saved (has to exist)
    % Return values:
    %  * hasChanged: true if file and saved hash do not match

    hash = getFileHash(sourceFilename);
    fid = fopen(hashFilename);
    tline = fgetl(fid);
    fclose(fid);
    hasChanged = ~strcmp(tline, hash);
end
