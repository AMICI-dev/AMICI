function [objectsstr, includesstr] = compileAMICIDependencies(wrap_path, objectFolder, o_suffix, COPT, DEBUG)
    %COMPILEAMICIDEPENDENCIES Compiles Sundials and SuiteSparse libraries required by AMICI

    sundials_path = fullfile(wrap_path,'sundials');
    sundials_ver = '2.7.0';
    
    ssparse_path = fullfile(wrap_path,'SuiteSparse');
    ssparse_ver = '4.5.3';
    
    lapack_path = fullfile(wrap_path,'lapack-3.5.0'); % currently not used, lapack implementation still needs to be done
    lapack_ver = '3.5.0';
    
    % compile directory
    if(~exist(objectFolder, 'dir'))
        mkdir(objectFolder)
    end
    
    version_file = fullfile(objectFolder, 'versions.txt');
    [del_sundials, del_ssparse, del_lapack] = checkVersions(version_file, sundials_ver, ssparse_ver, lapack_ver);
    
    % assemble objectsstr
    objectsstr = '';
    objects_sundials = getObjectsSundials(o_suffix);
    for j=1:length(objects_sundials)
        objectsstr = strcat(objectsstr,' "',fullfile(objectFolder,objects_sundials{j}),'"');
    end
    
    objects_ssparse = getObjectsSSparse(o_suffix);
    for j=1:length(objects_ssparse)
        objectsstr = strcat(objectsstr,' "',fullfile(objectFolder,objects_ssparse{j}),'"');
    end

    includesstr = getIncludeString(wrap_path, sundials_path, ssparse_path);
    
    % collect files that need to be recompiled
    sources_sundials = getSourcesSundials();
    sourcesToCompile = '';
    for j=1:length(sources_sundials)
        if(del_sundials || ~exist(fullfile(objectFolder,objects_sundials{j}), 'file'))
            sourcesToCompile = [sourcesToCompile, ' ', fullfile(sundials_path,sources_sundials{j})];
        end
    end
    sources_ssparse = getSourcesSSparse();
    for j=1:length(sources_ssparse)
        if(del_ssparse || ~exist(fullfile(objectFolder,objects_ssparse{j}), 'file'))
            sourcesToCompile = [sourcesToCompile, ' ', fullfile(ssparse_path,sources_ssparse{j})];
        end
    end
    
    % compile
    if(~strcmp(sourcesToCompile, ''))
        eval(['mex ' DEBUG ' ' COPT ' -c -outdir ' ...
            objectFolder ...
            includesstr ' ' sourcesToCompile ]);
    end
    
    % only write versions.txt if we are done compiling 
    fid = fopen(version_file,'w');
    fprintf(fid,[sundials_ver '\r']);
    fprintf(fid,[ssparse_ver '\r']);
    fprintf(fid,[lapack_ver '\r']);
    fclose(fid); 
end

function includesstr = getIncludeString(wrap_path, sundials_path, ssparse_path)
    includesstr = '';
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(wrap_path), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(wrap_path, 'src'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(wrap_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'KLU','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'AMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'COLAMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'BTF','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'SuiteSparse_config'), '"');
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
        objects_ssparse = strrep(objects_ssparse, '.o', o_suffix);
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
        objects_sundials = strrep(objects_sundials, '.o', o_suffix);
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
