function this = generateC(this)
    % generateC generates the c files which will be used in the compilation.
    %
    % Return values:
    %  this: model definition object @type amimodel

    if(strcmp(this.wtype,'iw'))
        xvec = 'x,';
        dxvec = 'dx,';
        sdxvec = 'sdx,';
        dxBvec = 'dxB,';
        rtcj = 'cj,';
    else
        xvec = '';
        dxvec = '';
        sdxvec = '';
        dxBvec = '';
        rtcj = '';
    end
    
    
    % write fun ccode
    
    for ifuns = 1:length(this.funs)
        if(this.fun.(this.funs{ifuns}))
            fprintf([this.funs{ifuns} ' | ']);
            fid = fopen(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' this.funs{ifuns} '.c']),'w');
            fprintf(fid,'\n');
            fprintf(fid,'#include <include/symbolic_functions.h>\n');
            fprintf(fid,'#include <string.h>\n');
            fprintf(fid,'#include <mex.h>\n');
            fprintf(fid,'#include <include/udata.h>\n');
            if(strcmp(this.funs{ifuns},'JBand'))
                fprintf(fid,['#include "' this.modelname '_J.h"\n']);
            elseif(strcmp(this.funs{ifuns},'JBandB'))
                fprintf(fid,['#include "' this.modelname '_JB.h"\n']);
            elseif(strcmp(this.funs{ifuns},'sxdot'))
                fprintf(fid,['#include "' this.modelname '_JSparse.h"\n']);
                fprintf(fid,['#include "' this.modelname '_dxdotdp.h"\n']);
            elseif(strcmp(this.funs{ifuns},'qBdot'))
                fprintf(fid,['#include "' this.modelname '_dxdotdp.h"\n']);
            end
            fprintf(fid,'\n');
            % function definition
            fprintf(fid,['int ' this.funs{ifuns} '_' this.modelname '' getFunArgs(this,this.funs{ifuns}) ' {\n']);
            if(strcmp(this.funs{ifuns},'JBand'))
                fprintf(fid,['return(J_' this.modelname '(N,t,' rtcj 'x,' dxvec 'xdot,J,user_data,tmp1,tmp2,tmp3));']);
            elseif(strcmp(this.funs{ifuns},'JBandB'))
                fprintf(fid,['return(JB_' this.modelname '(NeqBdot,t,' rtcj 'x,' dxvec 'xB,' dxBvec 'xdotB,J,user_data,tmp1,tmp2,tmp3));']);
            else
                fprintf(fid,'UserData udata = (UserData) user_data;\n');
                printlocalvars(this.funs{ifuns},this,fid);
                if( strcmp(this.funs{ifuns},'sx0') || strcmp(this.funs{ifuns},'sdx0') )
                    fprintf(fid,'switch (plist[ip]) {\n');
                    writeCcode_sensi(this, this.funs{ifuns}, fid);
                    fprintf(fid,'}\n');
                elseif( strcmp(this.funs{ifuns},'sxdot') )

                    if(strcmp(this.wtype,'iw'))
                        fprintf(fid,'int ip;\n');
                        fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                        fprintf(fid,'realtype *sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                        fprintf(fid,'realtype *sdx_tmp = N_VGetArrayPointer(sdx[plist[ip]]);\n');
                        fprintf(fid,'realtype *sxdot_tmp = N_VGetArrayPointer(sxdot[plist[ip]]);\n');
                        fprintf(fid,['memset(sxdot_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                        fprintf(fid,'switch (plist[ip]) {\n');
                        writeCcode_sensi(this, this.funs{ifuns}, fid);
                        fprintf(fid,'}\n');
                        fprintf(fid,'}\n');
                    else
                        fprintf(fid,'int status;\n');
                        fprintf(fid,'if(ip == 0) {\n');
                        fprintf(fid,['    status = JSparse_' this.modelname '(t,' rtcj 'x,' dxvec 'NULL,tmp_J,user_data,NULL,NULL,NULL);\n']);
                        fprintf(fid,['    status = dxdotdp_' this.modelname '(t,x,tmp_dxdotdp,user_data);\n']);
                        fprintf(fid,'}\n');
                        writeCcode(this, this.funs{ifuns}, fid);
                    end
                elseif( strcmp(this.funs{ifuns},'qBdot') )
                    fprintf(fid,'int status;\n');
                    fprintf(fid,'int ip;\n');
                    fprintf(fid,['status = dxdotdp_' this.modelname '(t,x,tmp_dxdotdp,user_data);\n']);
                    fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                    writeCcode(this, this.funs{ifuns}, fid);
                    fprintf(fid,'}\n');
                elseif( strcmp(this.funs{ifuns},'dydp')     || strcmp(this.funs{ifuns},'sy')         || ...
                        strcmp(this.funs{ifuns},'dtdp')       || strcmp(this.funs{ifuns},'drvaldp')  || strcmp(this.funs{ifuns},'sdeltadisc') || ...
                        strcmp(this.funs{ifuns},'dxdotdp')    || strcmp(this.funs{ifuns},'sroot')    || strcmp(this.funs{ifuns},'s2rootval')  || ...
                        strcmp(this.funs{ifuns},'s2root')     || strcmp(this.funs{ifuns},'srootval') || strcmp(this.funs{ifuns},'ideltadisc') || ...
                        strcmp(this.funs{ifuns},'dsigma_ydp') || strcmp(this.funs{ifuns},'dsigma_tdp')  )
                    fprintf(fid,'int ip;\n');
                    if(strcmp(this.funs{ifuns},'sdeltadisc'))
                        writeCcode(this, 'deltadisc', fid);
                    end
                    fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                    fprintf(fid,'switch (plist[ip]) {\n');
                    writeCcode_sensi(this, this.funs{ifuns}, fid);
                    fprintf(fid,'}\n');
                    fprintf(fid,'}\n');
                else
                    writeCcode(this, this.funs{ifuns}, fid);
                end
                printpostprocess(this.funs{ifuns},this,fid)
                fprintf(fid,'return(0);\n');
                fprintf(fid,'\n');
            end
            fprintf(fid,'}\n');
            fprintf(fid,'\n');
            fprintf(fid,'\n');
            fclose(fid);
        end
    end
    
    %
    %----------------------------------------------------------------
    % this.modelname.h
    % function declarations for this.modelname.c
    %----------------------------------------------------------------
    %
    
    fprintf('headers | ');
    fid = fopen(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '.h']),'w');
    fprintf(fid,['#ifndef _am_' this.modelname '_h\n']);
    fprintf(fid,['#define _am_' this.modelname '_h\n']);
    fprintf(fid,'\n');
    for ifuns = 1:length(this.funs)
        fprintf(fid,['#include "' this.modelname '_' this.funs{ifuns} '.h"\n']);
    end
    fprintf(fid,'\n');
    for ifuns = 1:length(this.funs)
        fprintf(fid,['int ' this.funs{ifuns} '_' this.modelname '' getFunArgs(this,this.funs{ifuns}) ';\n']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['#endif /* _am_' this.modelname '_h */\n']);
    fclose(fid);
    
    %
    %----------------------------------------------------------------
    % funname.h
    % function declarations for funname.c
    %----------------------------------------------------------------
    %
    
    for ifuns = 1:length(this.funs)
        fid = fopen(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' this.funs{ifuns} '.h']),'w');
        fprintf(fid,['#ifndef _am_' this.modelname '_' this.funs{ifuns} '_h\n']);
        fprintf(fid,['#define _am_' this.modelname '_' this.funs{ifuns} '_h\n']);
        fprintf(fid,'\n');
        fprintf(fid,['int ' this.funs{ifuns} '_' this.modelname '' getFunArgs(this,this.funs{ifuns}) ';\n']);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,['#endif /* _am_' this.modelname '_' this.funs{ifuns} '_h */\n']);
        fclose(fid);
    end
    %
    %----------------------------------------------------------------
    % wrapfunctions.c
    % function definitions
    %----------------------------------------------------------------
    %
    
    if(strcmp(this.wtype,'iw'))
        prefix = 'IDA';
        AMI = 'IDA';
        xvec = ' N_Vector x,';
        dxvec = ' N_Vector dx,';
        sdxvec = ' N_Vector sdx,';
        dxBvec = ' N_Vector dxB,';
        rtcj = ' realtype cj,';
        x = ', x';
        dx = ', dx';
        sdx = ', sdx';
        dxB = ', dxB';
        cj = ', cj';
        s = '*';
        one = '';
    else
        prefix = 'CV';
        AMI = 'CVode';
        xvec = '';
        dxvec = '';
        sdxvec = '';
        dxBvec = '';
        rtcj = '';
        x = '';
        dx = '';
        sdx = '';
        dxB = '';
        cj = '';
        s = '';
        one = '1';
    end
    
    
    fprintf('wrapfunctions | ');
    fid = fopen(fullfile(this.wrap_path,'models',this.modelname,'wrapfunctions.c'),'w');
    fprintf(fid,'                \n');
    fprintf(fid,'                \n');
    fprintf(fid,'                int wrap_init(void *cvode_mem, N_Vector x, N_Vector dx, realtype t){\n');
    fprintf(fid,['                    return ' AMI 'Init(cvode_mem, xdot_' this.modelname ', RCONST(t), x' dx ');\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_binit(void *cvode_mem, int which, N_Vector xB, N_Vector dxB, realtype t){\n');
    if(this.adjoint)
        fprintf(fid,['                    return ' AMI 'InitB(cvode_mem, which, xBdot_' this.modelname ', RCONST(t), xB' dxB ');\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_qbinit(void *cvode_mem, int which, N_Vector qBdot){\n');
    if(this.adjoint)
        fprintf(fid,['                    return ' AMI 'QuadInitB(cvode_mem, which, qBdot_' this.modelname ', qBdot);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SensInit1(void *cvode_mem, N_Vector *sx, N_Vector *sdx, void *user_data){\n');
    if(this.forward)
        fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
        fprintf(fid,['                    return ' AMI 'SensInit' one '(cvode_mem, np, sensi_meth, sxdot_' this.modelname ', sx' sdx ');\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                \n');
    fprintf(fid,'                int wrap_RootInit(void *cvode_mem, void *user_data){\n');
    fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
    fprintf(fid,['                    return ' AMI 'RootInit(cvode_mem, nr+ndisc, root_' this.modelname ');\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                \n');
    fprintf(fid,'                int fx0(N_Vector x0, void *user_data){\n');
    fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
    fprintf(fid,['                    return x0_' this.modelname '(x0, udata);\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdx0(N_Vector x0, N_Vector dx0, void *user_data){\n');
    if(strcmp(this.wtype,'iw'))
        fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
        fprintf(fid,['                    return dx0_' this.modelname '(x0, dx0, udata);\n']);
    else
        fprintf(fid,'                    return(0);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SetDenseJacFn(void *cvode_mem){\n');
    fprintf(fid,['                    return ' prefix 'DlsSetDenseJacFn(cvode_mem, J_' this.modelname ');\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SetSparseJacFn(void *cvode_mem){\n');
    fprintf(fid,['                    return ' prefix 'SlsSetSparseJacFn(cvode_mem, JSparse_' this.modelname ');\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SetBandJacFn(void *cvode_mem){\n');
    fprintf(fid,['                    return ' prefix 'DlsSetBandJacFn(cvode_mem, JBand_' this.modelname ');\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SetJacTimesVecFn(void *cvode_mem){\n');
    fprintf(fid,['                    return ' prefix 'SpilsSetJacTimesVecFn(cvode_mem, Jv_' this.modelname ');\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SetDenseJacFnB(void *cvode_mem,int which){\n');
    if(this.adjoint)
        fprintf(fid,['                    return ' prefix 'DlsSetDenseJacFnB(cvode_mem, which, JB_' this.modelname ');\n']);
    else
        fprintf(fid,'                    return(-1);\n');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SetSparseJacFnB(void *cvode_mem,int which){\n');
    if(this.adjoint)
        fprintf(fid,['                    return ' prefix 'SlsSetSparseJacFnB(cvode_mem, which, JSparseB_' this.modelname ');\n']);
    else
        fprintf(fid,'                    return(-1);\n');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SetBandJacFnB(void *cvode_mem,int which){\n');
    if(this.adjoint)
        fprintf(fid,['                    return ' prefix 'DlsSetBandJacFnB(cvode_mem, which, JBandB_' this.modelname ');\n']);
    else
        fprintf(fid,'                    return(-1);\n');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int wrap_SetJacTimesVecFnB(void *cvode_mem,int which){\n');
    if(this.adjoint)
        fprintf(fid,['                    return ' prefix 'SpilsSetJacTimesVecFnB(cvode_mem, which, JvB_' this.modelname ');\n']);
    else
        fprintf(fid,'                    return(-1);\n');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){\n');
    fprintf(fid,['                    return J_' this.modelname '(N,t' cj ',x' dx ',fx,J,user_data,tmp1,tmp2,tmp3);\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fJB(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){\n');
    if(this.adjoint)
        fprintf(fid,['                    return JB_' this.modelname '(N,t' cj ',x' dx ',xB' dxB ',fx,J,user_data,tmp1,tmp2,tmp3);\n']);
    else
        fprintf(fid,'                    return(-1);\n');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data){\n');
    fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
    fprintf(fid,'                    int ip;\n');
    fprintf(fid,'                    for (ip=0; ip<np; ip++) {\n');
    fprintf(fid,['                       sx0_' this.modelname '(plist[ip], sx0[ip]' x dx ', udata);\n']);
    fprintf(fid,'                    }\n');
    fprintf(fid,'                    return(0);\n');
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data){\n');
    if(strcmp(this.wtype,'iw'))
        fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
        fprintf(fid,'                    int ip;\n');
        fprintf(fid,'                    for (ip=0; ip<np; ip++) {\n');
        fprintf(fid,['                       sdx0_' this.modelname '(plist[ip], sdx0[ip]' x dx ', udata);\n']);
        fprintf(fid,'                    }\n');
        fprintf(fid,'                    return(0);\n');
    else
        fprintf(fid,'                    return(0);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data){\n');
    fprintf(fid,['                    return root_' this.modelname '(t, x' dx ', root, user_data);\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                \n');
    fprintf(fid,'                int fy(realtype t, int it, realtype *y, realtype *x, void *user_data){\n');
    fprintf(fid,['                    return y_' this.modelname '(t, it, y, x, user_data);\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'               int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data){\n');
    fprintf(fid,['                    return xdot_' this.modelname '(t,x' dx ',xdot,user_data);\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdxdotdp(realtype t, N_Vector x, realtype *dxdotdp, void *user_data){\n');
    if(~strcmp(this.wtype,'iw'))
        fprintf(fid,['                    return dxdotdp_' this.modelname '(t, x, dxdotdp, user_data);\n']);
    else
        fprintf(fid,'                    return(0);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdydp(realtype t, int it, realtype *dydp, realtype *y, realtype *x, void *user_data){\n');
    fprintf(fid,['                    return dydp_' this.modelname '(t, it, dydp, y, x, user_data);\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdydx(realtype t, int it, realtype *dydx, realtype *y, realtype *x, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return dydx_' this.modelname '(t, it, dydx, y, x, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdtdp(realtype t, realtype *dtdp, realtype *x, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return dtdp_' this.modelname '(t, dtdp, x, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdtdx(realtype t, realtype *dtdx, realtype *x, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return dtdx_' this.modelname '(t, dtdx, x, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdrvaldp(realtype t,realtype *drvaldp, realtype *x, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return drvaldp_' this.modelname '(t, drvaldp, x, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdrvaldx(realtype t,realtype *drvaldx, realtype *x, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return drvaldx_' this.modelname '(t, drvaldx, x, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fsy(realtype t, int it, realtype *sy, realtype *x, realtype *sx, void *user_data){\n');
    if(this.forward)
        fprintf(fid,['                    return sy_' this.modelname '(t, it, sy, x, sx, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fsroot(realtype t, int nroots, realtype *drvaldp, N_Vector x, N_Vector *sx, void *user_data){\n');
    if(this.forward)
        fprintf(fid,['                    return sroot_' this.modelname '(t, nroots, drvaldp, x, sx, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fsrootval(realtype t, int nroots, realtype *drootvaldp, N_Vector x, N_Vector *sx, void *user_data){\n');
    if(this.forward)
        fprintf(fid,['                    return srootval_' this.modelname '(t, nroots, drootvaldp, x, sx, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fs2root(realtype t, int nroots, realtype *ddrvaldpdp, N_Vector x, N_Vector *sx, void *user_data){\n');
    if(this.forward)
        fprintf(fid,['                    return s2root_' this.modelname '(t, nroots, ddrvaldpdp, x, sx, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fs2rootval(realtype t, int nroots, realtype *ddrootvaldpdp, N_Vector x, N_Vector *sx, void *user_data){\n');
    if(this.forward)
        fprintf(fid,['                    return s2rootval_' this.modelname '(t, nroots, ddrootvaldpdp, x, sx, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int deltadisc(realtype t, int idisc, N_Vector x, void *user_data){\n');
    if(~strcmp(this.wtype,'iw'))
        fprintf(fid,['                       return deltadisc_' this.modelname '(t, idisc, x, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    
    fprintf(fid,'                }\n');
    fprintf(fid,'                int bdeltadisc(realtype t, int idisc, N_Vector x, N_Vector xB, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                       return bdeltadisc_' this.modelname '(t, idisc, x, xB, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int ideltadisc(realtype t, int idisc, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                       return ideltadisc_' this.modelname '(t, idisc, x, xB, qBdot, user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int sdeltadisc(realtype t, int idisc, N_Vector x, N_Vector *sx, void *user_data){\n');
    if(this.forward)
        fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
        fprintf(fid,['                   return sdeltadisc_' this.modelname '(t, idisc, x, sx, udata);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fsigma_y(realtype t, realtype *sigma, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return sigma_y_' this.modelname '(t, sigma,user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdsigma_ydp(realtype t, realtype *dsigmadp, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return dsigma_ydp_' this.modelname '(t, dsigmadp,user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fsigma_t(realtype t, realtype *sigma, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return sigma_t_' this.modelname '(t, sigma,user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                int fdsigma_tdp(realtype t, realtype *dsigmadp, void *user_data){\n');
    if(this.adjoint)
        fprintf(fid,['                    return dsigma_tdp_' this.modelname '(t, dsigmadp,user_data);\n']);
    else
        fprintf(fid,'                    return(-1);');
    end
    fprintf(fid,'                }\n');
    fprintf(fid,'                \n');
    fprintf(fid,'                \n');
    fclose(fid);
    
    %
    %----------------------------------------------------------------
    % cvodewrapfunctions.h
    % function declarations for wrapfunctions.c
    %----------------------------------------------------------------
    %
    
    fid = fopen(fullfile(this.wrap_path,'models',this.modelname,'wrapfunctions.h'),'w');
    fprintf(fid,'#ifndef _am_wrapfunctions_h\n');
    fprintf(fid,'#define _am_wrapfunctions_h\n');
    fprintf(fid,'#include <math.h>\n');
    fprintf(fid,'#include <mex.h>\n');
    fprintf(fid,'\n');
    if(~strcmp(this.wtype,'iw'))
        fprintf(fid,'#include <include/cvodewrap.h>\n');
    else
        fprintf(fid,'#include <include/idawrap.h>\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,['#include "' this.modelname '.h"\n']);
    fprintf(fid,'\n');
    fprintf(fid,'#include <include/udata.h>\n');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'#define pi M_PI\n');
    fprintf(fid,'\n');
    fprintf(fid,'#include <wrapfunctions.c>\n');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'                int wrap_init(void *cvode_mem, N_Vector x, N_Vector dx, realtype t);\n');
    if(this.adjoint)
        fprintf(fid,'                int wrap_binit(void *cvode_mem, int which, N_Vector xB, N_Vector dxB, realtype t);\n');
        fprintf(fid,'                int wrap_qbinit(void *cvode_mem, int which, N_Vector qBdot);\n');
    end
    if(this.forward)
        fprintf(fid,'                int wrap_RootInit(void *cvode_mem, void *user_data);\n');
    end
    fprintf(fid,'                int fx0(N_Vector x0, void *user_data);\n');
    fprintf(fid,'                int fdx0(N_Vector x0, N_Vector dx0, void *user_data);\n');
    fprintf(fid,'                int wrap_SetDenseJacFn(void *cvode_mem);\n');
    fprintf(fid,'                int wrap_SetSparseJacFn(void *cvode_mem);\n');
    fprintf(fid,'                int wrap_SetBandJacFn(void *cvode_mem);\n');
    fprintf(fid,'                int wrap_SetJacTimesVecFn(void *cvode_mem);\n');
    fprintf(fid,'                int wrap_SetDenseJacFnB(void *cvode_mem,int which);\n');
    fprintf(fid,'                int wrap_SetSparseJacFnB(void *cvode_mem,int which);\n');
    fprintf(fid,'                int wrap_SetBandJacFnB(void *cvode_mem,int which);\n');
    fprintf(fid,'                int wrap_SetJacTimesVecFnB(void *cvode_mem,int which);\n');
    fprintf(fid,'                int fJ(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);\n');
    fprintf(fid,'                int fJB(long int N, realtype t, realtype cj, N_Vector x, N_Vector dx, N_Vector xB, N_Vector dxB, N_Vector fx, DlsMat J, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);\n');
    fprintf(fid,'                int fsx0(N_Vector *sx0, N_Vector x, N_Vector dx, void *user_data);\n');
    fprintf(fid,'                int fsdx0(N_Vector *sdx0, N_Vector x, N_Vector dx, void *user_data);\n');
    fprintf(fid,'                int froot(realtype t, N_Vector x, N_Vector dx, realtype *root, void *user_data);\n');
    fprintf(fid,'                int fy(realtype t, int it, realtype *y, realtype *x, void *user_data);\n');
    fprintf(fid,'                int fxdot(realtype t, N_Vector x, N_Vector dx, N_Vector xdot, void *user_data);\n');
    fprintf(fid,'                int fdxdotdp(realtype t, N_Vector x, realtype *dxdotdp, void *user_data);\n');
    fprintf(fid,'                int fdydp(realtype t, int it, realtype *dydp, realtype *y, realtype *x, void *user_data);\n');
    if(this.adjoint)
        fprintf(fid,'                int fdydx(realtype t, int it, realtype *dydx, realtype *y, realtype *x, void *user_data);\n');
    end
    fprintf(fid,'                int fdtdp(realtype t, realtype *dtdp, realtype *x, void *user_data);\n');
    fprintf(fid,'                int fdtdx(realtype t, realtype *drvaldp, realtype *x, void *user_data);\n');
    fprintf(fid,'                int fdrvaldp(realtype t, realtype *drvaldp, realtype *x, void *user_data);\n');
    fprintf(fid,'                int fdrvaldx(realtype t, realtype *drvaldx, realtype *x, void *user_data);\n');
    fprintf(fid,'                int fsigma_y(realtype t, realtype *sigma_y, void *user_data);\n');
    fprintf(fid,'                int fdsigma_ydp(realtype t, realtype *dsigma_ydp, void *user_data);\n');
    fprintf(fid,'                int fsigma_t(realtype t, realtype *sigma_t, void *user_data);\n');
    fprintf(fid,'                int fdsigma_tdp(realtype t, realtype *dsigma_tdp, void *user_data);\n');
    fprintf(fid,'#endif /* _LW_cvodewrapfunctions */\n');
    fclose(fid);
    
    fprintf('\r')
    
end

function printlocalvars(fun,struct,fid)
    % printlocalvars prints the C code for the initialisation of local variables into the file specified by fid.
    %
    % Parameters:
    %  fun: function name @type string
    %  struct: this struct must contain all necessary symbolic definitions @type symbolic
    %  fid: file id in which the final expression is written @type fileid
    %
    % Return values:
    %  Nothing
    
    
    nx = length(struct.sym.x);
    ny = length(struct.sym.y);
    nr = struct.nr;
    ndisc = struct.ndisc;
    
    switch(fun)
        case 'xdot'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'realtype *xdot_tmp = N_VGetArrayPointer(xdot);\n');
            fprintf(fid,['memset(xdot_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'xBdot'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'realtype *xB_tmp = N_VGetArrayPointer(xB);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dxB_tmp = N_VGetArrayPointer(dxB);\n');
            end
            fprintf(fid,'realtype *xBdot_tmp = N_VGetArrayPointer(xBdot);\n');
            fprintf(fid,['memset(xBdot_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'qBdot'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'realtype *xB_tmp = N_VGetArrayPointer(xB);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dxB_tmp = N_VGetArrayPointer(dxB);\n');
            end
            fprintf(fid,'realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);\n');
            fprintf(fid,'memset(qBdot_tmp,0,sizeof(realtype)*np);\n');
        case 'x0'
            fprintf(fid,'realtype *x0_tmp = N_VGetArrayPointer(x0);\n');
            fprintf(fid,['memset(x0_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'dx0'
            fprintf(fid,'realtype *x0_tmp = N_VGetArrayPointer(x0);\n');
            fprintf(fid,'realtype *dx0_tmp = N_VGetArrayPointer(dx0);\n');
            fprintf(fid,['memset(dx0_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'Jv'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'realtype *v_tmp = N_VGetArrayPointer(v);\n');
            fprintf(fid,'realtype *Jv_tmp = N_VGetArrayPointer(Jv);\n');
            fprintf(fid,['memset(Jv_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'JvB'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'realtype *xB_tmp = N_VGetArrayPointer(xB);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dxB_tmp = N_VGetArrayPointer(dxB);\n');
            end
            fprintf(fid,'realtype *vB_tmp = N_VGetArrayPointer(vB);\n');
            fprintf(fid,'realtype *JvB_tmp = N_VGetArrayPointer(JvB);\n');
            fprintf(fid,['memset(JvB_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'JBand'
            % nothing
        case 'J'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,['memset(J->data,0,sizeof(realtype)*' num2str(nx^2) ');\n']);
        case 'JSparse'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'SlsSetToZero(J);\n');
            for i = 1:length(struct.rowvals)
                fprintf(fid,['J->rowvals[' num2str(i-1) '] = ' num2str(struct.rowvals(i)) ';\n']);
            end
            for i = 1:length(struct.colptrs)
                fprintf(fid,['J->colptrs[' num2str(i-1) '] = ' num2str(struct.colptrs(i)) ';\n']);
            end
        case 'JBandB'
            % nothing
        case 'JB'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'realtype *xB_tmp = N_VGetArrayPointer(xB);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dxB_tmp = N_VGetArrayPointer(dxB);\n');
            end
            fprintf(fid,['  memset(JB->data,0,sizeof(realtype)*' num2str(nx^2) ');\n']);
        case 'JSparseB'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'  SlsSetToZero(JB);\n');
            for i = 1:length(struct.rowvalsB)
                fprintf(fid,['  JB->rowvals[' num2str(i-1) '] = ' num2str(struct.rowvalsB(i)) ';\n']);
            end
            for i = 1:length(struct.colptrsB)
                fprintf(fid,['  JB->colptrs[' num2str(i-1) '] = ' num2str(struct.colptrsB(i)) ';\n']);
            end
        case 'sxdot'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            if(~strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *sx_tmp = N_VGetArrayPointer(sx);\n');
                fprintf(fid,'realtype *sxdot_tmp = N_VGetArrayPointer(sxdot);\n');
                fprintf(fid,['memset(sxdot_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
            end
        case 'sx0'
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            end
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,'realtype *sx0_tmp = N_VGetArrayPointer(sx0);\n');
            fprintf(fid,['memset(sx0_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'sdx0'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            fprintf(fid,'realtype *sdx0_tmp = N_VGetArrayPointer(sdx0);\n');
            fprintf(fid,['memset(sdx0_tmp,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'y'
            
        case 'sy'
            
        case 'dydp'
            
        case 'dydx'
            
            fprintf(fid,['memset(dydx,0,sizeof(realtype)*' num2str(ny*nx) ');\n']);
        case 'dtdx'
            fprintf(fid,['memset(dtdx,0,sizeof(realtype)*' num2str(nr*nx) ');\n']);
        case 'drvaldx'
            fprintf(fid,['  memset(drvaldx,0,sizeof(realtype)*' num2str(nr*nx) ');\n']);
        case 'dtdp'
            
        case 'drvaldp'
            
        case 'root'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            if(strcmp(struct.wtype,'iw'))
                fprintf(fid,'realtype *dx_tmp = N_VGetArrayPointer(dx);\n');
            end
            fprintf(fid,['memset(root,0,sizeof(realtype)*' num2str(nr+ndisc) ');\n']);
        case 'sroot'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,'realtype *sx_tmp;\n');
            
        case 's2root'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,'realtype *sx_tmp;\n');
            
        case 'rootval'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,['memset(rootval,0,sizeof(realtype)*' num2str(nr+ndisc) ');\n']);
        case 'srootval'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,'realtype *sx_tmp;\n');
        case 's2rootval'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,'realtype *sx_tmp;\n');
        case 'deltadisc'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,['realtype deltadisc[' num2str(nx) '];\n']);
            fprintf(fid,['memset(deltadisc,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'bdeltadisc'
            fprintf(fid,'realtype *xB_tmp = N_VGetArrayPointer(xB);\n');
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,['realtype bdeltadisc[' num2str(nx) '];\n']);
            fprintf(fid,['memset(bdeltadisc,0,sizeof(realtype)*' num2str(nx) ');\n']);
        case 'ideltadisc'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,'realtype *xB_tmp = N_VGetArrayPointer(xB);\n');
            fprintf(fid,'realtype *qBdot_tmp = N_VGetArrayPointer(qBdot);\n');
            fprintf(fid,'realtype *ideltadisc;\n');
            fprintf(fid,'ideltadisc = mxMalloc(sizeof(realtype)*np);\n');
            fprintf(fid,'memset(ideltadisc,0,sizeof(realtype)*np);\n');
        case 'sdeltadisc'
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid,'realtype *sx_tmp;\n');
            fprintf(fid,['realtype deltadisc[' num2str(nx) '];\n']);
            fprintf(fid,'realtype *sdeltadisc;\n');
            fprintf(fid,['memset(deltadisc,0,sizeof(realtype)*' num2str(nx) ');\n']);
            fprintf(fid,['sdeltadisc = mxMalloc(sizeof(realtype)*' num2str(nx) '*np);\n']);
            fprintf(fid,['memset(sdeltadisc,0,sizeof(realtype)*' num2str(nx) '*np);\n']);
        case 'dxdotdp'
            fprintf(fid,['memset(dxdotdp,0,sizeof(realtype)*' num2str(nx) '*np);\n']);
            fprintf(fid,'realtype *x_tmp = N_VGetArrayPointer(x);\n');
        case 'sigma_y'
            fprintf(fid,['memset(sigma_y,0,sizeof(realtype)*' num2str(ny) ');\n']);
        case 'dsigma_ydp'
            fprintf(fid,['memset(dsigma_ydp,0,sizeof(realtype)*' num2str(ny) '*np);\n']);
        case 'sigma_t'
            fprintf(fid,['memset(sigma_t,0,sizeof(realtype)*' num2str(nr) ');\n']);
        case 'dsigma_tdp'
            fprintf(fid,['memset(dsigma_tdp,0,sizeof(realtype)*' num2str(nr) '*np);\n']);
        otherwise
            error(['unkown function: ' fun])
    end
    
end

function printpostprocess(fun,struct,fid)
    % printpostprocess prints the C code for the postprocessing into the file specified by fid.
    %
    % Parameters:
    %  fun: function name @type string
    %  struct: this struct must contain all necessary symbolic definitions @type symbolic
    %  fid: file id in which the final expression is written @type fileid
    %
    % Return values:
    %  Nothing
    
    nx = length(struct.sym.x);

    switch(fun)
        case 'xdot'
            %         ivar = 'ix';
            %         cvar = 'xdot_tmp';
            %         nvar = num2str(nx);
            %         printpp;
        case 'xBdot'
            %         ivar = 'ix';
            %         cvar = 'xBdot_tmp';
            %         nvar = num2str(nx);
            %         printpp;
        case 'qBdot'
            %         ivar = 'jp';
            %         cvar = 'qBdot_tmp';
            %         nvar = 'np';
            %         printpp;
        case 'x0'
            %         ivar = 'ix';
            %         cvar = 'x0_tmp';
            %         nvar =  num2str(nx);
            %         printpp;
        case 'dx0'
            %         ivar = 'ix';
            %         cvar = 'dx0_tmp';
            %         nvar =  num2str(nx);
            %         printpp;
        case 'Jv'
            %         ivar = 'ix';
            %         cvar = 'Jv_tmp';
            %         nvar =  num2str(nx);
            %         printpp;
        case 'JvB'
            %         ivar = 'ix';
            %         cvar = 'JvB_tmp';
            %         nvar =  num2str(nx);
            %         printpp;
        case 'JBand'
            % nothing
        case 'J'
            %         ivar = 'iJ';
            %         cvar = 'J->data';
            %         nvar =  num2str(nx^2);
            %         printpp;
        case 'JSparse'
            %         ivar = 'iJ';
            %         cvar = 'J->data';
            %         nvar =  num2str(length(struct.sparseidx));
            %         printpp;
        case 'JBBand'
            % nothing
        case 'JB'
            %         ivar = 'iJ';
            %         cvar = 'JB->data';
            %         nvar =  num2str(nx^2);
            %         printpp;
        case 'JSparseB'
            %         ivar = 'iJ';
            %         cvar = 'JB->data';
            %         nvar =  num2str(length(struct.sparseidxB));
            %         printpp;
        case 'sxdot'
            %         ivar = 'ix';
            %         cvar = 'sxdot_tmp';
            %         nvar =  num2str(nx);
            %         printpp;
        case 'sx0'
            %         ivar = 'ix';
            %         cvar = 'sx0_tmp';
            %         nvar =  num2str(nx);
            %         printpp;
        case 'sdx0'
            %         ivar = 'ix';
            %         cvar = 'sdx0_tmp';
            %         nvar =  num2str(nx);
            %         printpp;
        case 'y'
            %         ivar = 'iy';
            %         cvar = 'y';
            %         nvar =  num2str(ny);
            %         printpp;
        case 'sy'
            %
        case 'dydp'
            %
        case 'dydx'
            %
        case 'dtdx'
            %
        case 'drvaldx'
            %
        case 'dtdp'
            %
        case 'drvaldp'
            %
        case 'root'
            %
        case 'sroot'
            %
        case 's2root'
            %
        case 'rootval'
            %
        case 'srootval'
            %
        case 's2rootval'
            %
        case 'deltadisc'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['for(ix = 0; ix<' num2str(nx) ';ix++){;\n']);
            fprintf(fid,'  x_tmp[ix] += deltadisc[ix];\n');
            fprintf(fid,'};\n');
        case 'bdeltadisc'
            fprintf(fid,'int ix;\n');
            fprintf(fid,['for(ix = 0; ix<' num2str(nx) ';ix++){;\n']);
            fprintf(fid,'  xB_tmp[ix] += bdeltadisc[ix];\n');
            fprintf(fid,'};\n');
        case 'ideltadisc'
            fprintf(fid,'for(ip = 0; ip<np;ip++){;\n');
            fprintf(fid,'  qBdot_tmp[plist[ip]] += ideltadisc[plist[ip]];\n');
            fprintf(fid,'};\n');
            fprintf(fid,'mxFree(ideltadisc);\n');
        case 'sdeltadisc'
            fprintf(fid,'int ix;\n');
            fprintf(fid,'for(ip = 0; ip<np;ip++){\n');
            fprintf(fid,'  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
            fprintf(fid,['  for(ix = 0; ix<' num2str(nx) ';ix++){\n']);
            fprintf(fid,'    sx_tmp[ix] += sdeltadisc[plist[ip]+np*ix];\n');
            fprintf(fid,'  }\n');
            fprintf(fid,'}\n');
            fprintf(fid,['for(ix = 0; ix<' num2str(nx) ';ix++){\n']);
            fprintf(fid,'  x_tmp[ix] += deltadisc[ix];\n');
            fprintf(fid,'};\n');
            fprintf(fid,'mxFree(sdeltadisc);\n');
        case 'dxdotdp'
            % nothing
        case 'sigma_y'
            % nothing
        case 'dsigma_ydp'
            % nothing
        case 'sigma_t'
            % nothing
        case 'dsigma_tdp'
            % nothing
        otherwise
            error(['unkown function: ' fun])
    end
    
%     function printpp
%         fprintf(fid,['int ' ivar ';\n']);
%         fprintf(fid,['for (' ivar '=0; ' ivar '<' nvar '; ' ivar '++) {\n']);
%         fprintf(fid,['  if(mxIsNaN(' cvar '[' ivar '])) ' cvar '[' ivar '] = 0.0;\n']);
%         fprintf(fid,'}\n');
%     end
    
end