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

    [objectsstr, includesstr] = compileAMICIDependencies(this.wrap_path, o_suffix);
        
    includesstr = strcat(includesstr,' -I"', modelSourceFolder, '"');

    % append model object files
    for j=1:length(this.funs)
        objectsstr = strcat(objectsstr,...
            ' "',fullfile(modelSourceFolder, [this.modelname '_' this.funs{j} o_suffix]),'"');
    end

    % generate compile flags for the rest
    COPT = ['COPTIMFLAGS=''' this.coptim ' -DNDEBUG'' CXXFLAGS=''$CXXFLAGS -std=c++0x'''];
    if(this.debug)
        DEBUG = ' -g CXXFLAGS=''$CXXFLAGS -Wall  -std=c++0x -Wno-unused-function -Wno-unused-variable'' ';
        COPT = ''; % no optimization with debug flags!
    else
        DEBUG = '';
    end
    
    
    % generate hash for file and append debug string if we have an md5
    % file, check this hash against the contained hash
    cppsrc = {'symbolic_functions','spline','edata','rdata','udata','tdata', 'amici_interface_matlab', 'amici_misc'};
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
        if(this.recompile)
            recompile = 1;
        else
            recompile = checkHash(fullfile(modelSourceFolder,[this.modelname '_' this.funs{j}]),o_suffix,DEBUG);
        end
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
    
    if(this.cfun(1).(this.funs{j}))
        fprintf(['ffuns | ']);
        ffuns = cellfun(@(x) fullfile(modelSourceFolder,[this.modelname '_' x '.cpp']),this.funs,'UniformOutput',false);
        eval(['mex ' DEBUG COPT ...
            ' -c -outdir ' modelSourceFolder ' ' ...
            strrep(strcat(ffuns{:}),'.cpp','.cpp ') ' ' ...
            includesstr ...
            ' "' fullfile(amiciSourcePath,['symbolic_functions' o_suffix]) '"']);
        hash = getFileHash(fullfile(modelSourceFolder,[this.modelname '_' this.funs{j} '.cpp']));
        hash = [hash DEBUG];
        fid = fopen(...
            fullfile(modelSourceFolder,[this.modelname '_' this.funs{j} '_' mexext '.md5']...
            ),'w');
        fprintf(fid,hash);
        fclose(fid);
    end
    
    % compile the wrapfunctions object
    
    fprintf('wrapfunctions | '); 
    eval(['mex ' DEBUG COPT ...
        ' -c -outdir ' modelSourceFolder ' ' ...
        fullfile(modelSourceFolder,'wrapfunctions.cpp') ' ' ...
        includesstr]);
    
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
    
    % AMICI files which need to be recompiled, because of model-dependent
    % includes
    amiciSourceBaseNames = {'amici'};

    if(this.nxtrue == this.nx)
        amiciSourceList = '';
        for ii = 1:numel(amiciSourceBaseNames)
            amiciSourceList = [amiciSourceList, '"', fullfile(this.wrap_path, 'src', [amiciSourceBaseNames{ii}, '.cpp']), '" '];
        end
        eval(['mex ' DEBUG COPT ...
            ' -c -outdir ' fullfile(this.wrap_path,'src') ' ' ...
            amiciSourceList ' ' includesstr]);
    end
    
    amiciObjectList = '';
    for ii = 1:numel(amiciSourceBaseNames)
        amiciObjectList = [amiciObjectList, ' "', fullfile(this.wrap_path, 'src', [amiciSourceBaseNames{ii}, o_suffix]), '" '];
    end

    % Link object files
    eval(['mex ' DEBUG ' ' COPT ' ' CLIBS ...
        ' -output ' fullfile(modelSourceFolder,['ami_' this.modelname]) ...
        amiciObjectList ...
        ' "' fullfile(modelSourceFolder,['wrapfunctions' o_suffix]) '"' ...
        objectsstr ...
        includesstr ...
        ])
end    

function [del_sundials, del_ssparse, del_lapack] = checkVersions(version_file, sundials_ver, ssparse_ver, lapack_ver)
    % read version number from versions.txt and decide whether object have to
    % be regenerated
    if(~exist(version_file,'file'))
        del_sundials = true;
        del_ssparse = true;
        del_lapack = true;
        fid = fopen(version_file,'w');
        fprintf(fid,[sundials_ver '\r']);
        fprintf(fid,[ssparse_ver '\r']);
        fprintf(fid,[lapack_ver '\r']);
        fclose(fid);
    else
        fid = fopen(version_file,'r');
        sundials_objver = fgetl(fid);
        ssparse_objver = fgetl(fid);
        lapack_objver = fgetl(fid);
        fclose(fid);
        del_sundials = isnewer(sundials_ver,sundials_objver);
        if(del_sundials)
            display('Newer version of Sundials! Recompiling ...')
        end
        del_ssparse = isnewer(ssparse_ver,ssparse_objver);
        if(del_ssparse)
            display('Newer version of SuiteSparse! Recompiling ...')
        end
        del_lapack = isnewer(lapack_ver,lapack_objver);
        if(del_lapack)
            display('Newer version of Lapack! Recompiling ...')
        end
    end
end

function sources_sundials = getSourcesSundials()
    sources_sundials = {
        %    'src/cvodes/cvodes_lapack.c';
        fullfile('src','cvodes','cvodes_band.c');
        fullfile('src','cvodes','cvodes_bandpre.c');
        fullfile('src','cvodes','cvodes_bbdpre.c');
        fullfile('src','cvodes','cvodes_direct.c');
        fullfile('src','cvodes','cvodes_dense.c');
        fullfile('src','cvodes','cvodes_sparse.c');
        fullfile('src','cvodes','cvodes_diag.c');
        fullfile('src','cvodes','cvodea.c');
        fullfile('src','cvodes','cvodes.c');
        fullfile('src','cvodes','cvodes_io.c');
        fullfile('src','cvodes','cvodea_io.c');
        fullfile('src','cvodes','cvodes_spils.c');
        fullfile('src','cvodes','cvodes_spbcgs.c');
        fullfile('src','cvodes','cvodes_spgmr.c');
        fullfile('src','cvodes','cvodes_sptfqmr.c');
        fullfile('src','cvodes','cvodes_klu.c');
        fullfile('src','idas','idas.c');
        fullfile('src','idas','idas_sptfqmr.c');
        fullfile('src','idas','idas_spils.c');
        fullfile('src','idas','idas_spgmr.c');
        fullfile('src','idas','idas_spbcgs.c');
        fullfile('src','idas','idas_sparse.c');
        fullfile('src','idas','idas_klu.c');
        fullfile('src','idas','idas_io.c');
        fullfile('src','idas','idas_ic.c');
        fullfile('src','idas','idas_direct.c');
        fullfile('src','idas','idas_dense.c');
        fullfile('src','idas','idas_bbdpre.c');
        fullfile('src','idas','idas_band.c');
        fullfile('src','idas','idaa.c');
        fullfile('src','idas','idaa_io.c');
        fullfile('src','sundials','sundials_band.c');
        fullfile('src','sundials','sundials_dense.c');
        fullfile('src','sundials','sundials_sparse.c');
        fullfile('src','sundials','sundials_iterative.c');
        fullfile('src','sundials','sundials_nvector.c');
        fullfile('src','sundials','sundials_direct.c');
        fullfile('src','sundials','sundials_spbcgs.c');
        fullfile('src','sundials','sundials_spgmr.c');
        fullfile('src','sundials','sundials_sptfqmr.c');
        fullfile('src','sundials','sundials_math.c');
        fullfile('src','nvec_ser','nvector_serial.c');
        };
end

function sources_ssparse = getSourcesSSparse()
    sources_ssparse = {
        fullfile('KLU','Source','klu_analyze_given.c');
        fullfile('KLU','Source','klu_analyze.c');
        fullfile('KLU','Source','klu_defaults.c');
        fullfile('KLU','Source','klu_diagnostics.c');
        fullfile('KLU','Source','klu_dump.c');
        fullfile('KLU','Source','klu_extract.c');
        fullfile('KLU','Source','klu_factor.c');
        fullfile('KLU','Source','klu_free_numeric.c');
        fullfile('KLU','Source','klu_free_symbolic.c');
        fullfile('KLU','Source','klu_kernel.c');
        fullfile('KLU','Source','klu_memory.c');
        fullfile('KLU','Source','klu_refactor.c');
        fullfile('KLU','Source','klu_scale.c');
        fullfile('KLU','Source','klu_sort.c');
        fullfile('KLU','Source','klu_solve.c');
        fullfile('KLU','Source','klu_tsolve.c');
        fullfile('KLU','Source','klu.c');
        fullfile('AMD','Source','amd_1.c');
        fullfile('AMD','Source','amd_2.c');
        fullfile('AMD','Source','amd_aat.c');
        fullfile('AMD','Source','amd_control.c');
        fullfile('AMD','Source','amd_defaults.c');
        fullfile('AMD','Source','amd_dump.c');
        fullfile('AMD','Source','amd_global.c');
        fullfile('AMD','Source','amd_info.c');
        fullfile('AMD','Source','amd_order.c');
        fullfile('AMD','Source','amd_post_tree.c');
        fullfile('AMD','Source','amd_postorder.c');
        fullfile('AMD','Source','amd_preprocess.c');
        fullfile('AMD','Source','amd_valid.c');
        fullfile('COLAMD','Source','colamd.c');
        fullfile('BTF','Source','btf_maxtrans.c');
        fullfile('BTF','Source','btf_order.c');
        fullfile('BTF','Source','btf_strongcomp.c');
        fullfile('SuiteSparse_config','SuiteSparse_config.c');
        };
end

function objects_ssparse = getObjectsSSparse(o_suffix)
    
    objects_ssparse = {
        'klu_analyze_given.o';
        'klu_analyze.o';
        'klu_defaults.o';
        'klu_diagnostics.o';
        'klu_dump.o';
        'klu_extract.o';
        'klu_factor.o';
        'klu_free_numeric.o';
        'klu_free_symbolic.o';
        'klu_kernel.o';
        'klu_memory.o';
        'klu_refactor.o';
        'klu_scale.o';
        'klu_sort.o';
        'klu_solve.o';
        'klu_tsolve.o';
        'klu.o';
        'amd_1.o';
        'amd_2.o';
        'amd_aat.o';
        'amd_control.o';
        'amd_defaults.o';
        'amd_dump.o';
        'amd_global.o';
        'amd_info.o';
        'amd_order.o';
        'amd_post_tree.o';
        'amd_postorder.o';
        'amd_preprocess.o';
        'amd_valid.o';
        'colamd.o';
        'btf_maxtrans.o';
        'btf_order.o';
        'btf_strongcomp.o';
        'SuiteSparse_config.o';
        };
    
    if(~strcmp(o_suffix, '.o'))
        strrep(objects_ssparse, '.o', o_suffix);
    end
end

function objects_sundials = getObjectsSundials(o_suffix)
    objects_sundials = {
        %    'cvodes_lapack.o';
        'cvodes_band.o';
        'cvodes_bandpre.o';
        'cvodes_bbdpre.o';
        'cvodes_direct.o';
        'cvodes_dense.o';
        'cvodes_sparse.o';
        'cvodes_diag.o';
        'cvodea.o';
        'cvodes.o';
        'cvodes_io.o';
        'cvodea_io.o';
        'cvodes_spils.o';
        'cvodes_spbcgs.o';
        'cvodes_spgmr.o';
        'cvodes_sptfqmr.o';
        'cvodes_klu.o';
        'idas.o';
        'idas_sptfqmr.o';
        'idas_spils.o';
        'idas_spgmr.o';
        'idas_spbcgs.o';
        'idas_sparse.o';
        'idas_klu.o';
        'idas_io.o';
        'idas_ic.o';
        'idas_direct.o';
        'idas_dense.o';
        'idas_bbdpre.o';
        'idas_band.o';
        'idaa.o';
        'idaa_io.o';
        'sundials_band.o';
        'sundials_dense.o';
        'sundials_sparse.o';
        'sundials_iterative.o';
        'sundials_nvector.o';
        'sundials_direct.o';
        'sundials_spbcgs.o';
        'sundials_spgmr.o';
        'sundials_sptfqmr.o';
        'sundials_math.o';
        'nvector_serial.o';
        };
    
    if(~strcmp(o_suffix, '.o'))
        strrep(objects_sundials, '.o', o_suffix);
    end
end

function result = isnewer(ver1str,ver2str)
    % isnewer checks whether the version indicated in ver1str is newer than
    % the on in ver2str
    %
    % Parameters:
    %  ver1str: version string 1, this should be a string in the format %d.%d.%d @type string
    %  ver2str: version string 1, this should be a string in the format %d.%d.%d @type string
    %
    % Return values:
    %  result: flag indicating whether ver1 ist newer than ver2 @type boolean
    
    ver1Parts = getParts(ver1str);
    ver2Parts = getParts(ver2str);
    if ver2Parts(1) ~= ver1Parts(1)     % major version
        result = ver2Parts(1) < ver1Parts(1);
    elseif ver2Parts(2) ~= ver1Parts(2) % minor version
        result = ver2Parts(2) < ver1Parts(2);
    elseif ver2Parts(3) ~= ver1Parts(3) % revision version
        result = ver2Parts(3) < ver1Parts(3);
    else
        result = ver2Parts(4) < ver1Parts(4);
    end
end 

function parts = getParts(V)
    % getParts takes an input version string and returns an array
    % containing the major minor and subminor version number
    %
    % Parameters:
    %  V: version string, this should be a string in the format %d.%d.%d @type string
    %
    % Return values:
    %  parts: array containing the version numbers @type double
    
    parts = sscanf(V, '%d.%d.%d.%d')';
    if length(parts) < 3
        parts(3) = 0; % zero-fills to 3 elements
    end
    if length(parts) < 4
        parts(4) = 0; % zero-fills to 3 elements
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


function [objectsstr, includesstr] = compileAMICIDependencies(wrap_path, o_suffix)

    sundials_path = fullfile(wrap_path,'sundials');
    sundials_ver = '2.7.0';
    
    ssparse_path = fullfile(wrap_path,'SuiteSparse');
    ssparse_ver = '4.5.3';
    
    lapack_path = fullfile(wrap_path,'lapack-3.5.0'); % currently not used, lapack implementation still needs to be done
    lapack_ver = '3.5.0';
    
    % compile directory
    if(~exist(fullfile(wrap_path,'models',mexext), 'dir'))
        mkdir(fullfile(wrap_path,'models',mexext))
    end
    
    version_file = fullfile(wrap_path,'models',mexext,'versions.txt');
    [del_sundials, del_ssparse, del_lapack] = checkVersions(version_file, sundials_ver, ssparse_ver, lapack_ver);
    
    % assemble objectsstr
    objectsstr = '';
    objects_sundials = getObjectsSundials(o_suffix);
    for j=1:length(objects_sundials)
        objectsstr = strcat(objectsstr,' "',fullfile(wrap_path,'models',mexext,objects_sundials{j}),'"');
    end
    
    objects_ssparse = getObjectsSSparse(o_suffix);
    for j=1:length(objects_ssparse)
        objectsstr = strcat(objectsstr,' "',fullfile(wrap_path,'models',mexext,objects_ssparse{j}),'"');
    end
        
    includesstr = '';
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'src','cvodes'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(wrap_path), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(wrap_path, 'src'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(wrap_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'KLU','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'AMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'COLAMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'BTF','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'SuiteSparse_config'), '"');
    
    % compile all the sundials objects if we haven't done so yet
    sources_sundials = getSourcesSundials();
    for j=1:length(sources_sundials)
        if(~exist(fullfile(wrap_path,'models',mexext,objects_sundials{j}), 'file') || del_sundials)
            eval(['mex COPTIMFLAGS=''-O3 -DNDEBUG'' -c -outdir '...
                fullfile(wrap_path,'models',mexext) ...
                includesstr ' ' ...
                fullfile(sundials_path,sources_sundials{j})]);
        end
    end
    
    % compile all the suitesparse objects if we haven't done so yet
    sources_ssparse = getSourcesSSparse();
    for j=1:length(sources_ssparse)
        if(~exist(fullfile(wrap_path,'models',mexext,objects_ssparse{j}), 'file') || del_ssparse)
            eval(['mex COPTIMFLAGS=''-O3 -DNDEBUG'' -c -outdir ' ...
                fullfile(wrap_path,'models',mexext) ...
                includesstr ' ' ...
                fullfile(ssparse_path,sources_ssparse{j})]);
        end
    end
    
    % only write versions.txt if we are done compiling 
    fid = fopen(fullfile(wrap_path,'models',mexext,'versions.txt'),'w');
    fprintf(fid,[sundials_ver '\r']);
    fprintf(fid,[ssparse_ver '\r']);
    fprintf(fid,[lapack_ver '\r']);
    fclose(fid);
end
