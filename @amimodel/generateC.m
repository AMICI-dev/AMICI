function generateC(this)
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
    
    for ifun = this.funs
        if(isfield(this.fun,ifun{1}))
            fprintf([ifun{1} ' | ']);
            fid = fopen(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' ifun{1} '.c']),'w');
            fprintf(fid,'\n');
            fprintf(fid,'#include <include/symbolic_functions.h>\n');
            fprintf(fid,'#include <string.h>\n');
            fprintf(fid,'#include <mex.h>\n');
            if( strfind(this.fun.(ifun{1}).argstr,'user_data') )
                fprintf(fid,'#include <include/udata.h>\n');
            end
            if( strfind(this.fun.(ifun{1}).argstr,'temp_data') )
                fprintf(fid,'#include <include/tdata.h>\n');
                % avoid conflicts
                fprintf(fid,'#undef t\n');
                fprintf(fid,'#undef x\n');
                fprintf(fid,'#undef x_tmp\n');
                fprintf(fid,'#undef dzdp\n'); 
                fprintf(fid,'#undef dzdx\n'); 
                fprintf(fid,'#undef dsigma_zdp\n'); 
            end
            if(strcmp(ifun{1},'JBand'))
                fprintf(fid,['#include "' this.modelname '_J.h"\n']);
            elseif(strcmp(ifun{1},'JBandB'))
                fprintf(fid,['#include "' this.modelname '_JB.h"\n']);
            elseif(strcmp(ifun{1},'sxdot'))
                fprintf(fid,['#include "' this.modelname '_JSparse.h"\n']);
                fprintf(fid,['#include "' this.modelname '_dxdotdp.h"\n']);
            elseif(strcmp(ifun{1},'qBdot'))
                fprintf(fid,['#include "' this.modelname '_dxdotdp.h"\n']);
            end
            fprintf(fid,'\n');
            % function definition
            fprintf(fid,['int ' ifun{1} '_' this.modelname '' this.fun.(ifun{1}).argstr ' {\n']);
            if(strcmp(ifun{1},'JBand'))
                fprintf(fid,['return(J_' this.modelname removeTypes(this.fun.J.argstr) ');']);
            elseif(strcmp(ifun{1},'JBandB'))
                fprintf(fid,['return(JB_' this.modelname removeTypes(this.fun.JB.argstr) ');']);
            else
                if( strfind(this.fun.(ifun{1}).argstr,'user_data') )
                    fprintf(fid,'UserData udata = (UserData) user_data;\n');
                end
                if( strfind(this.fun.(ifun{1}).argstr,'temp_data') )
                    fprintf(fid,'TempData tdata = (TempData) temp_data;\n');
                end
                this.fun.(ifun{1}).printLocalVars(this,fid);
                if( strcmp(ifun{1},'sxdot') )
                    if(strcmp(this.wtype,'iw'))
                        fprintf(fid,'int ip;\n');
                        fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                        fprintf(fid,'realtype *sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                        fprintf(fid,'realtype *sdx_tmp = N_VGetArrayPointer(sdx[plist[ip]]);\n');
                        fprintf(fid,'realtype *sxdot_tmp = N_VGetArrayPointer(sxdot[plist[ip]]);\n');
                        fprintf(fid,['memset(sxdot_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                        fprintf(fid,'switch (plist[ip]) {\n');
                        this.fun.(ifun{1}).writeCcode_sensi(this,fid);
                        fprintf(fid,'}\n');
                        fprintf(fid,'}\n');
                    else
                        fprintf(fid,'int status;\n');
                        fprintf(fid,'if(ip == 0) {\n');
                        fprintf(fid,['    status = JSparse_' this.modelname '(t,' rtcj 'x,' dxvec 'xdot,tmp_J,user_data,NULL,NULL,NULL);\n']);
                        fprintf(fid,['    status = dxdotdp_' this.modelname '(t,tmp_dxdotdp,x,user_data);\n']);
                        fprintf(fid,'}\n');
                        this.fun.(ifun{1}).writeCcode(this,fid);
                    end
                elseif( strcmp(ifun{1},'qBdot') )
                    fprintf(fid,'int status;\n');
                    fprintf(fid,'int ip;\n');
                    fprintf(fid,['status = dxdotdp_' this.modelname '(t,tmp_dxdotdp,x,user_data);\n']);
                    fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                    this.fun.(ifun{1}).writeCcode(this,fid);
                    fprintf(fid,'}\n');
                elseif(this.fun.(ifun{1}).sensiflag)
                    fprintf(fid,'int ip;\n');
                    fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                    if(ismember('*sx',this.fun.(ifun{1}).nvecs))
                        fprintf(fid,'sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                    end
                    if(ismember('*sx0',this.fun.(ifun{1}).nvecs))
                        fprintf(fid,'sx0_tmp = N_VGetArrayPointer(sx0[plist[ip]]);\n');
                        fprintf(fid,['memset(sx0_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                    end
                    fprintf(fid,'switch (plist[ip]) {\n');
                    this.fun.(ifun{1}).writeCcode_sensi(this,fid);
                    fprintf(fid,'}\n');
                    fprintf(fid,'}\n');
                else
                    this.fun.(ifun{1}).writeCcode(this,fid);
                end
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
    for ifun = this.funs
        fprintf(fid,['#include "' this.modelname '_' ifun{1} '.h"\n']);
    end
    fprintf(fid,'\n');
    for ifun = this.funs
        if(isfield(this.fun,ifun{1}))
            fun = this.fun.(ifun{1});
        else
            fun = amifun(ifun{1},this);
        end
        fprintf(fid,['int ' ifun{1} '_' this.modelname '' fun.argstr ';\n']);
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
    
    for ifun = this.funs
        if(isfield(this.fun,ifun{1}))
            fun = this.fun.(ifun{1});
        else
            fun = amifun(ifun{1},this);
        end
        fid = fopen(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' ifun{1} '.h']),'w');
        fprintf(fid,['#ifndef _am_' this.modelname '_' ifun{1} '_h\n']);
        fprintf(fid,['#define _am_' this.modelname '_' ifun{1} '_h\n']);
        fprintf(fid,'\n');
        fprintf(fid,['int ' ifun{1} '_' this.modelname '' fun.argstr ';\n']);
        fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,['#endif /* _am_' this.modelname '_' ifun{1} '_h */\n']);
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
    fprintf(fid,'#include "wrapfunctions.h"\n');
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
    fprintf(fid,['                    return ' AMI 'RootInit(cvode_mem, ne, root_' this.modelname ');\n']);
    fprintf(fid,'                }\n');
    fprintf(fid,'                \n');
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
    
    ffuns = {'x0','dx0','sx0','sdx0','J','JB','root',...
        'y','sy','dydp','dydx','z','sz','sz_tf','dzdp','dzdx',...
        'xdot','xBdot','qBdot','dxdotdp','deltax','deltasx','deltaxB','deltaqB',...
        'sigma_y','dsigma_ydp','sigma_z','dsigma_zdp',...
        'Jy','Jz','dJydx','dJydp','dJzdx','dJzdp','sJy','sJz'};
    
    for iffun = ffuns
        % check whether the function was generated, otherwise generate (but
        % whithout symbolic expressions)
        if(isfield(this.fun,iffun{1}))
            fun = this.fun.(iffun{1});
        else
            fun = amifun(iffun{1},this);
        end
        fprintf(fid,['                int f' iffun{1} fun.fargstr '{\n']);
        % if the function was generated, we can return it, otherwise return
        % an error
        if(ismember(iffun{1},this.funs))
            fprintf(fid,['                    return ' iffun{1} '_' this.modelname removeTypes(fun.argstr) ';\n']);
        else
            if(strcmp(iffun{1},'sx0') || strcmp(iffun{1},'dx0') || strcmp(iffun{1},'sdx0'))
                % these two should always work, if they are not required
                % they should act as placeholders
                fprintf(fid,'                    return(0);\n');
            else
                fprintf(fid,'                    return(-1);\n');
            end
        end
        fprintf(fid,'                }\n');
        fprintf(fid,'                \n');
    end
    
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
    fprintf(fid,'\n');
    fprintf(fid,'                int wrap_init(void *cvode_mem, N_Vector x, N_Vector dx, realtype t);\n');
    fprintf(fid,'                int wrap_binit(void *cvode_mem, int which, N_Vector xB, N_Vector dxB, realtype t);\n');
    fprintf(fid,'                int wrap_qbinit(void *cvode_mem, int which, N_Vector qBdot);\n');
    fprintf(fid,'                int wrap_RootInit(void *cvode_mem, void *user_data);\n');
    fprintf(fid,'                int wrap_SensInit1(void *cvode_mem, N_Vector *sx, N_Vector *sdx, void *user_data);\n');
    fprintf(fid,'                int wrap_SetDenseJacFn(void *cvode_mem);\n');
    fprintf(fid,'                int wrap_SetSparseJacFn(void *cvode_mem);\n');
    fprintf(fid,'                int wrap_SetBandJacFn(void *cvode_mem);\n');
    fprintf(fid,'                int wrap_SetJacTimesVecFn(void *cvode_mem);\n');
    fprintf(fid,'                int wrap_SetDenseJacFnB(void *cvode_mem,int which);\n');
    fprintf(fid,'                int wrap_SetSparseJacFnB(void *cvode_mem,int which);\n');
    fprintf(fid,'                int wrap_SetBandJacFnB(void *cvode_mem,int which);\n');
    fprintf(fid,'                int wrap_SetJacTimesVecFnB(void *cvode_mem,int which);\n');
    for iffun = ffuns
        % check whether the function was generated, otherwise generate (but
        % whithout symbolic expressions)
        if(isfield(this.fun,iffun{1}))
            fun = this.fun.(iffun{1});
        else
            fun = amifun(iffun{1},this);
        end
        fprintf(fid,['                int f' iffun{1} fun.fargstr ';\n']);
    end
    fprintf(fid,'#endif /* _LW_cvodewrapfunctions */\n');
    fclose(fid);
    
    fprintf('\r')
end

    
function argstr = removeTypes(argstr)
    % removeTypes transforms an argument string from a string with variable
    % types (for function definition) to a string without variable types
    % (for function calling)
    %
    % Parameters:
    %  argstr: function definition argument string @type *char
    %
    % Return values:
    %  argstr: function call argument string @type *char
    
    argstr = strrep(argstr,'realtype','');
    argstr = strrep(argstr,'int','');
    argstr = strrep(argstr,'void','');
    argstr = strrep(argstr,'N_Vector','');
    argstr = strrep(argstr,'long','');
    argstr = strrep(argstr,'DlsMat','');
    argstr = strrep(argstr,'SlsMat','');
    argstr = strrep(argstr,'*','');
    argstr = strrep(argstr,' ','');
    argstr = strrep(argstr,',',', ');
    
end