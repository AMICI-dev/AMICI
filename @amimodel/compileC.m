function compileC(this)
    % compileC compiles the mex simulation file
    %
    % Return values:
    %  this: model definition object @type amimodel
    
    if(ispc)
        o_suffix = '.obj';
    else
        o_suffix = '.o';
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
    
    [objectsstr, includesstr] = compileAMICIDependencies(this.wrap_path, o_suffix, COPT, DEBUG);
        
    includesstr = strcat(includesstr,' -I"', modelSourceFolder, '"');

    % append model object files
    for j=1:length(this.funs)
        objectsstr = strcat(objectsstr,...
            ' "',fullfile(modelSourceFolder, [this.modelname '_' this.funs{j} o_suffix]),'"');
    end    
    
    % generate hash for file and append debug string if we have an md5
    % file, check this hash against the contained hash
    cppsrc = {'amici', 'symbolic_functions','spline','edata','rdata','udata','tdata', 'amici_interface_matlab', 'amici_misc', 'amici_model_functions'};
    for srcfile = cppsrc
        baseFilename = fullfile(amiciSourcePath,srcfile{1});
        
        recompile = this.recompile || checkHash(baseFilename,o_suffix,DEBUG);
       
        if(recompile)
            fprintf([srcfile{1} ' | ']);
            eval(['mex ' DEBUG COPT ...
                ' -c -outdir ' amiciSourcePath ...
                includesstr ' '...
                [baseFilename '.cpp']]);
            updateFileHash(baseFilename, DEBUG);
        end
        objectsstr = strcat(objectsstr,' "', [baseFilename o_suffix],'"');
    end  
    
    % do the same for all the this.funs
    for j=1:length(this.funs)
        baseFilename = fullfile(modelSourceFolder,[this.modelname '_' this.funs{j}]);

        recompile = this.recompile || checkHash(baseFilename,o_suffix,DEBUG);

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
    
    % compile the wrapfunctions object
    
    fprintf('wrapfunctions | '); 
    eval(['mex ' DEBUG COPT ...
        ' -c -outdir ' modelSourceFolder ' ' ...
        fullfile(modelSourceFolder,'wrapfunctions.cpp') ' ' ...
        includesstr]);
    objectsstr = [objectsstr, ' "' fullfile(modelSourceFolder,['wrapfunctions' o_suffix]) '"'];

    % now we have compiled everything model specific, so we can replace hashes.mat to prevent recompilation
    try
    movefile(fullfile(modelSourceFolder,'hashes_new.mat'),...
        fullfile(modelSourceFolder,'hashes.mat'),'f');
    end
    

    fprintf('amici | ');

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

    % Link object files
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
