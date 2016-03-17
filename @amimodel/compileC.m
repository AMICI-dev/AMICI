function compileC(this)
    % compileC compiles the mex simulation file
    %

    
    sundials_path = fullfile(this.wrap_path,'sundials-2.6.2');
    sundials_ver = '2.6.2.1';
    
    ssparse_path = fullfile(this.wrap_path,'SuiteSparse');
    ssparse_ver = '4.4.4';
    
    lapack_path = fullfile(this.wrap_path,'lapack-3.5.0'); % currently not used, lapack implementation still needs to be done
    lapack_ver = '3.5.0';
    
    includesstr = '';
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'src','cvodes'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(this.wrap_path, 'models', this.modelname ), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(this.wrap_path), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(this.wrap_path, 'src'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(this.wrap_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'KLU','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'AMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'COLAMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'BTF','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'SuiteSparse_config'), '"');
    
    % compile directory
    if(~exist(fullfile(this.wrap_path,'models',mexext), 'dir'))
        mkdir(fullfile(this.wrap_path,'models',mexext))
    end
    
    % read version number from versions.txt and decide whether object have to
    % be regenerated
    if(~exist(fullfile(this.wrap_path,'models',mexext,'versions.txt'),'file'))
        del_sundials = true;
        del_ssparse = true;
        del_lapack = true;
        fid = fopen(fullfile(this.wrap_path,'models',mexext,'versions.txt'),'w');
        fprintf(fid,[sundials_ver '\r']);
        fprintf(fid,[ssparse_ver '\r']);
        fprintf(fid,[lapack_ver '\r']);
        fclose(fid);
    else
        fid = fopen(fullfile(this.wrap_path,'models',mexext,'versions.txt'),'r');
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
        fid = fopen(fullfile(this.wrap_path,'models',mexext,'versions.txt'),'w');
        fprintf(fid,[sundials_ver '\r']);
        fprintf(fid,[ssparse_ver '\r']);
        fprintf(fid,[lapack_ver '\r']);
        fclose(fid);
    end
    
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
    sourcesstr_sundials = '';
    for j=1:length(sources_sundials)
        sourcesstr_sundials = strcat(sourcesstr_sundials,' "',fullfile(sundials_path,sources_sundials{j}),'"');
    end
    
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
        fullfile('COLAMD','Source','colamd_global.c');
        fullfile('COLAMD','Source','colamd.c');
        fullfile('BTF','Source','btf_maxtrans.c');
        fullfile('BTF','Source','btf_order.c');
        fullfile('BTF','Source','btf_strongcomp.c');
        fullfile('SuiteSparse_config','SuiteSparse_config.c');
        };
    sourcesstr_ssparse = '';
    for j=1:length(sources_ssparse)
        sourcesstr_ssparse = strcat(sourcesstr_ssparse,' "', fullfile(ssparse_path,sources_ssparse{j}),'"');
    end
    
    sourcesstr_funs = '';
    for j=1:length(this.funs)
        sourcesstr_funs = strcat(sourcesstr_funs,...
            ' "', fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' this.funs{j} '.c']),'"');
    end
    
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
    if(ispc)
        objects_sundials = strrep(objects_sundials, '.o', '.obj');
    end
    
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
        'colamd_global.o';
        'colamd.o';
        'btf_maxtrans.o';
        'btf_order.o';
        'btf_strongcomp.o';
        'SuiteSparse_config.o';
        };
    if(ispc)
        objects_ssparse = strrep(objects_ssparse, '.o', '.obj');
    end
    
    
    % assemble objectsstr
    objectsstr = '';
    for j=1:length(objects_sundials)
        objectsstr = strcat(objectsstr,' "',fullfile(this.wrap_path,'models',mexext,objects_sundials{j}),'"');
    end
    
    for j=1:length(objects_ssparse)
        objectsstr = strcat(objectsstr,' "',fullfile(this.wrap_path,'models',mexext,objects_ssparse{j}),'"');
    end
    
    for j=1:length(this.funs)
        if(ispc)
            o_suffix = '.obj';
        else
            o_suffix = '.o';
        end
        objectsstr = strcat(objectsstr,...
            ' "',fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' this.funs{j} o_suffix]),'"');
    end
    
    objectsstr = strcat(objectsstr,' "',fullfile(this.wrap_path,'src',['symbolic_functions' o_suffix]),'"');
    objectsstr = strcat(objectsstr,' "',fullfile(this.wrap_path,'src',['amici' o_suffix]),'"');
    
    
    % compile all the sundials objects if we haven't done so yet
    for j=1:length(sources_sundials)
        if(~exist(fullfile(this.wrap_path,'models',mexext,objects_sundials{j}), 'file') || del_sundials)
            eval(['mex COPTIMFLAGS=''-O3 -DNDEBUG'' -c -outdir '...
                fullfile(this.wrap_path,'models',mexext) ...
                includesstr ' ' ...
                fullfile(sundials_path,sources_sundials{j})]);
        end
    end
    
    % compile all the sundials objects if we haven't done so yet
    for j=1:length(sources_ssparse)
        if(~exist(fullfile(this.wrap_path,'models',mexext,objects_ssparse{j}), 'file') || del_ssparse)
            eval(['mex COPTIMFLAGS=''-O3 -DNDEBUG'' -c -outdir ' ...
                fullfile(this.wrap_path,'models',mexext) ...
                includesstr ' ' ...
                fullfile(ssparse_path,sources_ssparse{j})]);
        end
    end
      
    
    % generate compile flags for the rest
    COPT = ['COPTIMFLAGS=''' this.coptim ' -DNDEBUG'''];
    if(this.debug)
        DEBUG = '-g';
        COPT = ''; % no optimization with debug flags!
    else
        DEBUG = '';
    end
    
    
    % generate hash for file and append debug string if we have an md5
    % file, check this hash against the contained hash
    if(this.recompile)
        recompile = 1;
    else
        recompile = checkHash(fullfile(this.wrap_path,'src','symbolic_functions'),o_suffix,DEBUG);
    end
    if(recompile)
        fprintf('symbolic_functions | ');
        eval(['mex ' DEBUG COPT ...
            ' -c -outdir ' fullfile(this.wrap_path,'src') ...
            includesstr ' ' ...
            fullfile(this.wrap_path,'src','symbolic_functions.c')]);
        hash = getFileHash(fullfile(this.wrap_path,'src','symbolic_functions.c'));
        hash = [hash DEBUG];
        fid = fopen(fullfile(this.wrap_path,'src',['symbolic_functions' '_' mexext '.md5']),'w');
        fprintf(fid,hash);
        fclose(fid);
    end
    
    
    
    
    % do the same for all the this.funs
    for j=1:length(this.funs)
        if(this.recompile)
            recompile = 1;
        else
            recompile = checkHash(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' this.funs{j}]),o_suffix,DEBUG);
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
    if(this.cfun(1).dxdotdp)
        this.cfun(1).sxdot = 1;
        this.cfun(1).qBdot = 1;
    end
    if(this.cfun(1).w)
        this.recompile = 1;
    end
    
    recompileWrapFunction = false; 
    % if any of the functions in this.funs is recompiled, we also need to
    % recompile the wrapfunction object
     
    for j=1:length(this.funs)
        if(this.cfun(1).(this.funs{j}))
            recompileWrapFunction = true;
            fprintf([this.funs{j} ' | ']);
            eval(['mex ' DEBUG COPT ...
                ' -c -outdir ' fullfile(this.wrap_path,'models',this.modelname) ...
                includesstr ...
                ' "' fullfile(this.wrap_path,'src',['symbolic_functions' o_suffix]) '" ' ...
                fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' this.funs{j} '.c'])]);
            hash = getFileHash(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' this.funs{j} '.c']));
            hash = [hash DEBUG];
            fid = fopen(...
                fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' this.funs{j} '_' mexext '.md5']...
                ),'w');
            fprintf(fid,hash);
            fclose(fid);
        end
    end
    
    % compile the wrapfunctions object
    
    if(recompileWrapFunction || ~exist(fullfile(this.wrap_path,'models',this.modelname,['wrapfunctions' o_suffix]),'file'))
        fprintf('wrapfunctions | ');
        eval(['mex ' DEBUG COPT ...
                ' -c -outdir ' fullfile(this.wrap_path,'models',this.modelname) ...
                includesstr ' '...
                fullfile(this.wrap_path,'models',this.modelname,'wrapfunctions.c')]);
    end
    

    fprintf('amici | ');
    eval(['mex ' DEBUG COPT ...
        ' -c -outdir ' fullfile(this.wrap_path,'src') ...
        includesstr ' ' ...
        fullfile(this.wrap_path,'src','amici.c')]);
    
    if(isunix)
        if(~ismac)
            CLIBS = 'CLIBS=''\$CLIBS -lrt''';
        else
            CLIBS = [];
        end
    else
        CLIBS = [];
    end
    
    if(this.nxtrue ~= this.nx)
        cstr = 'amiwrapo2.c';
    else
        cstr = 'amiwrap.c';
    end

    prefix = 'ami';
    
    eval(['mex ' DEBUG ' ' COPT ' ' CLIBS ...
        ' -output ' fullfile(this.wrap_path,'models',this.modelname,[prefix '_' this.modelname]) ...
        ' ' fullfile(this.wrap_path,cstr)  ...
        includesstr ...
        objectsstr ...
        ' "',fullfile(this.wrap_path,'models',this.modelname,['wrapfunctions' o_suffix]),'"'
        ])

    
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
    else                                  % revision version
        result = ver2Parts(3) < ver1Parts(3);
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
    
    parts = sscanf(V, '%d.%d.%d')';
    if length(parts) < 3
        parts(3) = 0; % zero-fills to 3 elements
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
    
    
    
function recompile = checkHash(filestr,o_suffix,DEBUG)
    % checkHash checks whether filestr.c  has already been compiled as
    % filestr.o and whether the md5 hash of filestr.c matches the one in
    % filestr.md5
    %
    % Parameters:
    %  filestr: path of the file @type string
    %  o_suffix: OS specific suffix for compiled objects
    %  DEBUG: debug flag
    %
    % Return values:
    %  recompile: flag indicating whether we need to recompile filestr.c
    
    if(~exist([filestr o_suffix],'file'))
        % file does not exist, we need to recompile
        recompile = 1;
    else
        if(~exist([filestr '_' mexext '.md5'], 'file'))
            % hash does not exist, we need to recompile
            recompile = 1;
        else
            hash = getFileHash([filestr '.c']);
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