function amiActivateDebugMode( flag )
    %amiActivateDebugMode recompiles sundials objects in order to
    %enable/disable debugging of sundials specific functions
    
    this.wrap_path=fileparts(fileparts(mfilename('fullpath')));  
    sundials_path = fullfile(this.wrap_path,'sundials-2.6.2');
    
    ssparse_path = fullfile(this.wrap_path,'SuiteSparse');
    ssparse_ver = '4.4.4';
    
    includesstr = '';
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'src','cvodes'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'src','idas'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(this.wrap_path), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(this.wrap_path, 'src'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(this.wrap_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'KLU','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'AMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'COLAMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'BTF','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'SuiteSparse_config'), '"');
    
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
    
    if(flag)
        DEBUG = '-g';
    else
        DEBUG = '';
    end
    
    for j=1:length(sources_sundials)
        eval(['mex ' DEBUG ' COPTIMFLAGS=''-O3 -DNDEBUG'' -c -outdir '...
            fullfile(this.wrap_path,'models',mexext) ...
            includesstr ' ' ...
            fullfile(sundials_path,sources_sundials{j})]);
    end
    
    
end

