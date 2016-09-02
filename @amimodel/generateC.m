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
    dxvec = 'NULL,';
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
            fprintf(fid,'#undef dx\n');
            fprintf(fid,'#undef sigma_y\n');
            fprintf(fid,'#undef sigma_z\n');
            fprintf(fid,'#undef dsigma_ydp\n');
            fprintf(fid,'#undef dsigma_zdp\n');
        end
        if(strcmp(ifun{1},'JBand'))
            fprintf(fid,['#include "' this.modelname '_J.h"\n']);
        elseif(strcmp(ifun{1},'JBandB'))
            fprintf(fid,['#include "' this.modelname '_JB.h"\n']);
        elseif(strcmp(ifun{1},'sxdot'))
            fprintf(fid,['#include "' this.modelname '_JSparse.h"\n']);
            fprintf(fid,['#include "' this.modelname '_dxdotdp.h"\n']);
            if(strcmp(this.wtype,'iw'))
                fprintf(fid,['#include "' this.modelname '_dfdx.h"\n']);
                fprintf(fid,['#include "' this.modelname '_M.h"\n']);
            end
        elseif(strcmp(ifun{1},'dxdotdp') || strcmp(ifun{1},'qBdot'))
            fprintf(fid,['#include "' this.modelname '_dwdp.h"\n']);
        elseif( strcmp(ifun{1},'J') || strcmp(ifun{1},'JB') || strcmp(ifun{1},'JSparse') || strcmp(ifun{1},'JSparseB') || strcmp(ifun{1},'xBdot') || strcmp(ifun{1},'dfdx'))
            fprintf(fid,['#include "' this.modelname '_dwdx.h"\n']);
        elseif(strcmp(ifun{1},'qBdot'))
            fprintf(fid,['#include "' this.modelname '_dxdotdp.h"\n']);
        end
        fprintf(fid,['#include "' this.modelname '_w.h"\n']);
        fprintf(fid,'\n');
        % function definition
        fprintf(fid,['int ' ifun{1} '_' this.modelname '' this.fun.(ifun{1}).argstr ' {\n']);
        fprintf(fid,'int status = 0;\n');
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
            if(~isempty(strfind(this.fun.(ifun{1}).argstr,'N_Vector x')) && ~isempty(strfind(this.fun.(ifun{1}).argstr,'realtype t')))
                if(or(not(strcmp(this.wtype,'iw')),~isempty(strfind(this.fun.(ifun{1}).argstr,'N_Vector dx'))))
                    if(~any(ismember(ifun{1},{'w','sxdot','dxdotdp','dfdx','qBdot'})))
                        fprintf(fid,['status = w_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                    end
                end
            end
            if( strcmp(ifun{1},'sxdot') )
                if(strcmp(this.wtype,'iw'))
                    fprintf(fid,'int ip;\n');
                    fprintf(fid,['status = dfdx_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                    fprintf(fid,['status = M_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                    fprintf(fid,['status = dxdotdp_' this.modelname '(t,tmp_dxdotdp,x,' dxvec 'user_data);\n']);
                    fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                    fprintf(fid,'sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                    fprintf(fid,'sdx_tmp = N_VGetArrayPointer(sdx[plist[ip]]);\n');
                    fprintf(fid,'sxdot_tmp = N_VGetArrayPointer(sxdot[plist[ip]]);\n');
                    fprintf(fid,['memset(sxdot_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                    this.fun.(ifun{1}).writeCcode(this,fid);
                    fprintf(fid,'}\n');
                else
                    fprintf(fid,'if(ip == 0) {\n');
                    fprintf(fid,['    status = JSparse_' this.modelname '(t,' rtcj 'x,xdot,tmp_J,user_data,NULL,NULL,NULL);\n']);
                    fprintf(fid,['    status = dxdotdp_' this.modelname '(t,tmp_dxdotdp,x,' dxvec 'user_data);\n']);
                    fprintf(fid,'}\n');
                    this.fun.(ifun{1}).writeCcode(this,fid);
                end
            elseif( strcmp(ifun{1},'qBdot') )
                fprintf(fid,['status = dwdp_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                fprintf(fid,'switch (plist[ip]) {\n');
                this.fun.(ifun{1}).writeCcode_sensi(this,fid);
                fprintf(fid,'}\n');
                fprintf(fid,'}\n');
            elseif(this.fun.(ifun{1}).sensiflag)
                if( strcmp(ifun{1},'dxdotdp'))
                    fprintf(fid,['status = dwdp_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                end
                fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                if(ismember('*sx',this.fun.(ifun{1}).nvecs))
                    fprintf(fid,'sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                end
                if(ismember('*sx0',this.fun.(ifun{1}).nvecs))
                    fprintf(fid,'sx0_tmp = N_VGetArrayPointer(sx0[plist[ip]]);\n');
                    fprintf(fid,['memset(sx0_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                end
                if(ismember('*sdx0',this.fun.(ifun{1}).nvecs))
                    fprintf(fid,'sdx0_tmp = N_VGetArrayPointer(sdx0[plist[ip]]);\n');
                    fprintf(fid,['memset(sdx0_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                end
                fprintf(fid,'switch (plist[ip]) {\n');
                this.fun.(ifun{1}).writeCcode_sensi(this,fid);
                fprintf(fid,'}\n');
                fprintf(fid,'}\n');
            else
                if( strcmp(ifun{1},'J') || strcmp(ifun{1},'JB') || strcmp(ifun{1},'JSparse') || strcmp(ifun{1},'JSparseB') || strcmp(ifun{1},'xBdot') || strcmp(ifun{1},'dfdx') )
                    fprintf(fid,['status = dwdx_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                end
                this.fun.(ifun{1}).writeCcode(this,fid);
            end
            if(strcmp(ifun{1},'dxdotdp'))
                fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                fprintf(fid,['   for(ix = 0; ix<' num2str(this.nx) '; ix++) {\n']);
                fprintf(fid,['       if(mxIsNaN(dxdotdp[ix+ip*' num2str(this.nx) '])) {\n']);
                fprintf(fid,['           dxdotdp[ix+ip*' num2str(this.nx) '] = 0;\n']);
                fprintf(fid,'           if(!udata->am_nan_dxdotdp) {\n');
                fprintf(fid,'               mexWarnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'               udata->am_nan_dxdotdp = TRUE;\n');
                fprintf(fid,'           }\n');
                fprintf(fid,'       }\n');
                fprintf(fid,['       if(mxIsInf(dxdotdp[ix+ip*' num2str(this.nx) '])) {\n']);
                fprintf(fid,'           mexWarnMsgIdAndTxt("AMICI:mex:fdxdotdp:Inf","AMICI encountered an Inf value in dxdotdp, aborting.");\n');
                fprintf(fid,'           return(-1);\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'}\n');
            end
            if(this.nx>0)
                if(strcmp(ifun{1},'xdot'))
                    fprintf(fid,['for(ix = 0; ix<' num2str(this.nx) '; ix++) {\n']);
                    fprintf(fid,'   if(mxIsNaN(xdot_tmp[ix])) {\n');
                    fprintf(fid,'       xdot_tmp[ix] = 0;\n');
                    fprintf(fid,'       if(!udata->am_nan_xdot) {\n');
                    fprintf(fid,'           mexWarnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                    fprintf(fid,'           udata->am_nan_xdot = TRUE;\n');
                    fprintf(fid,'       }\n');
                    fprintf(fid,'   }\n');
                    fprintf(fid,'   if(mxIsInf(xdot_tmp[ix])) {\n');
                    fprintf(fid,'       mexWarnMsgIdAndTxt("AMICI:mex:fxdot:Inf","AMICI encountered an Inf value in xdot! Aborting simulation ... ");\n');
                    fprintf(fid,'       return(-1);\n');
                    fprintf(fid,'   }');
                    fprintf(fid,'   if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) {\n');
                    fprintf(fid,'       xdot_tmp[ix] = -xdot_tmp[ix];\n');
                    fprintf(fid,'   }\n');
                    fprintf(fid,'}\n');
                end
            end
            if(strcmp(ifun{1},'JSparse'))
                fprintf(fid,['for(inz = 0; inz<' num2str(this.nnz) '; inz++) {\n']);
                fprintf(fid,'   if(mxIsNaN(J->data[inz])) {\n');
                fprintf(fid,'       J->data[inz] = 0;\n');
                fprintf(fid,'       if(!udata->am_nan_JSparse) {\n');
                fprintf(fid,'           mexWarnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           udata->am_nan_JSparse = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'   if(mxIsInf(J->data[inz])) {\n');
                fprintf(fid,'       mexWarnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'J'))
                fprintf(fid,['for(ix = 0; ix<' num2str(this.nx^2) '; ix++) {\n']);
                fprintf(fid,'   if(mxIsNaN(J->data[ix])) {\n');
                fprintf(fid,'       J->data[ix] = 0;\n');
                fprintf(fid,'       if(!udata->am_nan_J) {\n');
                fprintf(fid,'           mexWarnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           udata->am_nan_J = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'   if(mxIsInf(J->data[ix])) {\n');
                fprintf(fid,'       mexWarnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'xBdot'))
                fprintf(fid,['for(ix = 0; ix<' num2str(this.nx) '; ix++) {\n']);
                fprintf(fid,'   if(mxIsNaN(xBdot_tmp[ix])) {\n');
                fprintf(fid,'       xBdot_tmp[ix] = 0;');
                fprintf(fid,'       if(!udata->am_nan_xBdot) {\n');
                fprintf(fid,'           mexWarnMsgIdAndTxt("AMICI:mex:fxBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           udata->am_nan_xBdot = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }');
                fprintf(fid,'   if(mxIsInf(xBdot_tmp[ix])) {\n');
                fprintf(fid,'       mexWarnMsgIdAndTxt("AMICI:mex:fxBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'qBdot'))
                fprintf(fid,'for(ip = 0; ip<np*ng; ip++) {\n');
                fprintf(fid,'   if(mxIsNaN(qBdot_tmp[ip])) {\n');
                fprintf(fid,'       qBdot_tmp[ip] = 0;');
                fprintf(fid,'       if(!udata->am_nan_qBdot) {\n');
                fprintf(fid,'           mexWarnMsgIdAndTxt("AMICI:mex:fqBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           udata->am_nan_qBdot = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }');
                fprintf(fid,'   if(mxIsInf(qBdot_tmp[ip])) {\n');
                fprintf(fid,'       mexWarnMsgIdAndTxt("AMICI:mex:fqBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }');
                fprintf(fid,'}\n');
            end
            fprintf(fid,'return(status);\n');
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
fprintf(fid,'                void init_modeldims(void *user_data){\n');
fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
fprintf(fid,['                   nx = ' num2str(this.nx) ';\n']);
fprintf(fid,['                   nxtrue = ' num2str(this.nxtrue) ';\n']);
fprintf(fid,['                   ny = ' num2str(this.ny) ';\n']);
fprintf(fid,['                   nytrue = ' num2str(this.nytrue) ';\n']);
fprintf(fid,['                   nz = ' num2str(this.nz) ';\n']);
fprintf(fid,['                   nztrue = ' num2str(this.nztrue) ';\n']);
fprintf(fid,['                   ne = ' num2str(this.nevent) ';\n']);
fprintf(fid,['                   ng = ' num2str(this.ng) ';\n']);
fprintf(fid,['                   nw = ' num2str(this.nw) ';\n']);
fprintf(fid,['                   ndwdx = ' num2str(this.ndwdx) ';\n']);
fprintf(fid,['                   ndwdp = ' num2str(this.ndwdp) ';\n']);
fprintf(fid,['                   nnz = ' num2str(this.nnz) ';\n']);
fprintf(fid,['                   ubw = ' num2str(this.ubw) ';\n']);
fprintf(fid,['                   lbw = ' num2str(this.lbw) ';\n']);
fprintf(fid,'                }\n');
fprintf(fid,'                int wrap_init(void *cvode_mem, N_Vector x, N_Vector dx, realtype t){\n');
fprintf(fid,['                    return ' AMI 'Init(cvode_mem, xdot_' this.modelname ', RCONST(t), x' dx ');\n']);
fprintf(fid,'                }\n');
fprintf(fid,'                int wrap_binit(void *cvode_mem, int which, N_Vector xB, N_Vector dxB, realtype t){\n');
if(this.adjoint)
    fprintf(fid,['                    return ' AMI 'InitB(cvode_mem, which, xBdot_' this.modelname ', RCONST(t), xB' dxB ');\n']);
else
    fprintf(fid,'                    return(-1);\n');
end
fprintf(fid,'                }\n');
fprintf(fid,'                int wrap_qbinit(void *cvode_mem, int which, N_Vector qBdot){\n');
if(this.adjoint)
    fprintf(fid,['                    return ' AMI 'QuadInitB(cvode_mem, which, qBdot_' this.modelname ', qBdot);\n']);
else
    fprintf(fid,'                    return(-1);\n');
end
fprintf(fid,'                }\n');
fprintf(fid,'                int wrap_SensInit1(void *cvode_mem, N_Vector *sx, N_Vector *sdx, void *user_data){\n');
if(this.forward)
    fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
    fprintf(fid,['                    return ' AMI 'SensInit' one '(cvode_mem, np, sensi_meth, sxdot_' this.modelname ', sx' sdx ');\n']);
else
    fprintf(fid,'                    return(-1);\n');
end
fprintf(fid,'                }\n');
fprintf(fid,'                \n');
fprintf(fid,'                int wrap_RootInit(void *cvode_mem, void *user_data){\n');
fprintf(fid,'                    UserData udata = (UserData) user_data;\n');
fprintf(fid,['                    return ' AMI 'RootInit(cvode_mem, ' num2str(this.nevent) ', root_' this.modelname ');\n']);
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

ffuns = {'x0','dx0','sx0','sdx0','J','JB','root','sroot','s2root','stau',...
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
% wrapfunctions.h
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
fprintf(fid,'                void init_modeldims(void *user_data);\n');
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