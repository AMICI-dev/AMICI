function generateCStandalone(this)
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
                fprintf(fid,['       if(amiIsNaN(dxdotdp[ix+ip*' num2str(this.nx) '])) {\n']);
                fprintf(fid,['           dxdotdp[ix+ip*' num2str(this.nx) '] = 0;\n']);
                fprintf(fid,'           if(!udata->am_nan_dxdotdp) {\n');
                fprintf(fid,'               warnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'               udata->am_nan_dxdotdp = TRUE;\n');
                fprintf(fid,'           }\n');
                fprintf(fid,'       }\n');
                fprintf(fid,['       if(amiIsInf(dxdotdp[ix+ip*' num2str(this.nx) '])) {\n']);
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fdxdotdp:Inf","AMICI encountered an Inf value in dxdotdp, aborting.");\n');
                fprintf(fid,'           return(-1);\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'}\n');
            end
            if(this.nx>0)
                if(strcmp(ifun{1},'xdot'))
                    fprintf(fid,['for(ix = 0; ix<' num2str(this.nx) '; ix++) {\n']);
                    fprintf(fid,'   if(amiIsNaN(xdot_tmp[ix])) {\n');
                    fprintf(fid,'       xdot_tmp[ix] = 0;\n');
                    fprintf(fid,'       if(!udata->am_nan_xdot) {\n');
                    fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                    fprintf(fid,'           udata->am_nan_xdot = TRUE;\n');
                    fprintf(fid,'       }\n');
                    fprintf(fid,'   }\n');
                    fprintf(fid,'   if(amiIsInf(xdot_tmp[ix])) {\n');
                    fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fxdot:Inf","AMICI encountered an Inf value in xdot! Aborting simulation ... ");\n');
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
                fprintf(fid,'   if(amiIsNaN(J->data[inz])) {\n');
                fprintf(fid,'       J->data[inz] = 0;\n');
                fprintf(fid,'       if(!udata->am_nan_JSparse) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           udata->am_nan_JSparse = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'   if(amiIsInf(J->data[inz])) {\n');
                fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'J'))
                fprintf(fid,['for(ix = 0; ix<' num2str(this.nx^2) '; ix++) {\n']);
                fprintf(fid,'   if(amiIsNaN(J->data[ix])) {\n');
                fprintf(fid,'       J->data[ix] = 0;\n');
                fprintf(fid,'       if(!udata->am_nan_J) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           udata->am_nan_J = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'   if(amiIsInf(J->data[ix])) {\n');
                fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'xBdot'))
                fprintf(fid,['for(ix = 0; ix<' num2str(this.nx) '; ix++) {\n']);
                fprintf(fid,'   if(amiIsNaN(xBdot_tmp[ix])) {\n');
                fprintf(fid,'       xBdot_tmp[ix] = 0;');
                fprintf(fid,'       if(!udata->am_nan_xBdot) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fxBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           udata->am_nan_xBdot = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }');
                fprintf(fid,'   if(amiIsInf(xBdot_tmp[ix])) {\n');
                fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fxBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'qBdot'))
                fprintf(fid,'for(ip = 0; ip<np; ip++) {\n');
                fprintf(fid,'   if(amiIsNaN(qBdot_tmp[ip])) {\n');
                fprintf(fid,'       qBdot_tmp[ip] = 0;');
                fprintf(fid,'       if(!udata->am_nan_qBdot) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fqBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           udata->am_nan_qBdot = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }');
                fprintf(fid,'   if(amiIsInf(qBdot_tmp[ip])) {\n');
                fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fqBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");\n');
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
fprintf(fid,['                   ny = ' num2str(this.ny) ';\n']);
fprintf(fid,['                   nytrue = ' num2str(this.nytrue) ';\n']);
fprintf(fid,['                   nz = ' num2str(this.nz) ';\n']);
fprintf(fid,['                   nztrue = ' num2str(this.nztrue) ';\n']);
fprintf(fid,['                   ne = ' num2str(this.nevent) ';\n']);
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

fprintf('CMakeLists | ');
generateCMakeFile(this);

fprintf('main | ');
generateMainC(this);

fprintf('\r');
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


function generateCMakeFile(this)
    CMakeFileName = fullfile(this.wrap_path,'models',this.modelname,'CMakeLists.txt');
    fid = fopen(CMakeFileName,'w');
    fprintf(fid, 'project(%s)\n', this.modelname);
    fprintf(fid, 'cmake_minimum_required(VERSION 2.8)\n');
    fprintf(fid, 'set(cmake_build_type Debug)\n');
    fprintf(fid, 'set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-unused-function")\n');
    fprintf(fid, 'add_definitions(-DAMICI_WITHOUT_MATLAB)\n');
    
    % sources
    fprintf(fid, '\nset(SRC_LIST\n');
    for f = {'main.c', 'wrapfunctions.c', ...
            fullfile(this.wrap_path, 'src/symbolic_functions.c'), ...
            fullfile(this.wrap_path, 'src/amici.c')
            }
        fprintf(fid, '%s\n', f{1});
    end
    for j=1:length(this.funs)
        %if(this.cfun(1).(this.funs{j}))
             funcName = this.funs{j};
             fprintf(fid, '%s_%s.c\n', this.modelname, funcName);
       % end
    end
    fprintf(fid, ')\n\n');
    
    %includes
    sundials_path = fullfile(this.wrap_path,'sundials-2.6.2');
    ssparse_path = fullfile(this.wrap_path,'SuiteSparse');
    
    includeDirs = {fullfile(this.wrap_path, 'models', this.modelname), ...
        fullfile(this.wrap_path),  ...
        fullfile(sundials_path, 'include'),  ...
        fullfile(ssparse_path, 'KLU', 'Include'), ...
        fullfile(ssparse_path, 'AMD', 'Include'), ...
        fullfile(ssparse_path, 'SuiteSparse_config'), ...
        fullfile(ssparse_path, 'COLAMD', 'Include'), ...
        fullfile(ssparse_path, 'BTF', 'Include')
    };
    for d = includeDirs
        fprintf(fid, 'include_directories("%s")\n', d{1});
    end
    
    fprintf(fid, '\n');
    fprintf(fid, 'add_executable(${PROJECT_NAME} ${SRC_LIST})\n\n');

    %libraries
    fprintf(fid, 'target_link_libraries(${PROJECT_NAME}\n');
    sundials_lib_path = '/home/dweindl/src/_libs/sundials-2.6.2/build/';
    ssparse_lib_path = '/home/dweindl/src/_libs/SuiteSparse/';
    libs = {
        fullfile(sundials_lib_path, 'install/lib/libsundials_nvecserial.so'), ...
        fullfile(sundials_lib_path, 'src/cvodes/libsundials_cvodes.so'), ...
        fullfile(ssparse_lib_path, 'lib/libcolamd.so'), ...
        fullfile(ssparse_lib_path, 'KLU/Lib/libklu.a'), ...
        fullfile(ssparse_lib_path, 'BTF/Lib/libbtf.a'), ...
        fullfile(ssparse_lib_path, 'AMD/Lib/libamd.a'), ...
        fullfile(ssparse_lib_path, 'COLAMD/Lib/libcolamd.a'), ...
        fullfile(ssparse_lib_path, 'SuiteSparse_config/libsuitesparseconfig.a'), ...
        '-lm'
    };
    for l = libs
        fprintf(fid, '"%s"\n', l{1});
    end
    fprintf(fid, ')\n');
    fclose(fid);
end

function generateMainC(this)
    mainFileName = fullfile(this.wrap_path,'models',this.modelname,'main.c');
    fid = fopen(mainFileName,'w');
    
    fprintf(fid, '#include <stdio.h>\n');
    fprintf(fid, '#include <stdlib.h>\n');
    fprintf(fid, '#include <string.h>\n');
    fprintf(fid, '#include "wrapfunctions.h" /* user functions */\n');
    fprintf(fid, '#include "include/symbolic_functions.h"\n');
    fprintf(fid, '#include <include/amici.h> /* amici functions */\n');
    fprintf(fid, '\n');
    fprintf(fid, 'UserData getUserData();\n');
    fprintf(fid, 'void processUserData(UserData udata);\n');
    fprintf(fid, 'ExpData getExperimentalData(UserData udata);\n');
    fprintf(fid, 'ReturnData getSimulationResults(UserData udata, ExpData edata);\n');
    fprintf(fid, 'ReturnData initReturnData(void *user_data, double *pstatus);\n');
    fprintf(fid, 'void freeTempDataAmiMem(UserData udata, TempData tdata, void *ami_mem, booleantype setupBdone, double *pstatus, int status);\n');
    fprintf(fid, 'void processReturnData(ReturnData rdata, UserData udata);\n');
    fprintf(fid, 'void printReturnData(ReturnData rdata, UserData udata);\n');
    fprintf(fid, 'void freeExpData(ExpData edata);\n');
    fprintf(fid, 'void freeReturnData(ReturnData rdata);\n');
    fprintf(fid, 'void freeUserData(UserData udata);\n');
    fprintf(fid, '\n');
    fprintf(fid, 'int main(int argc, char **argv)\n');
    fprintf(fid, '{  \n');
    fprintf(fid, '    UserData udata = getUserData();\n');
    fprintf(fid, '    if (udata == NULL) {\n');
    fprintf(fid, '        return 1;\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    ExpData edata = getExperimentalData(udata);\n');
    fprintf(fid, '    if (edata == NULL) {\n');
    fprintf(fid, '        freeUserData(udata);\n');
    fprintf(fid, '        return 1;\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    ReturnData rdata = getSimulationResults(udata, edata);\n');
    fprintf(fid, '    if (rdata == NULL) {\n');
    fprintf(fid, '        freeExpData(edata);\n');
    fprintf(fid, '        freeUserData(udata);\n');
    fprintf(fid, '        return 1;\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    processReturnData(rdata, udata);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    freeExpData(edata);\n');
    fprintf(fid, '    freeUserData(udata);\n');
    fprintf(fid, '    freeReturnData(rdata);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    return 0;\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'UserData getUserData() { // was: udata = setupUserData(prhs);\n');
    fprintf(fid, '    UserData udata; /* User udata structure */\n');
    fprintf(fid, '    udata = (UserData) malloc(sizeof *udata);\n');
    fprintf(fid, '    if (udata == NULL)\n');
    fprintf(fid, '        return(NULL);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    init_modeldims(udata);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* BEGIN SET USER DATA */\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* time vector, matlab: first argument */\n');
    fprintf(fid, '    nt = 1; // TODO\n');
    fprintf(fid, '    ts = linSpaceAlloc(0, 1, nt); // TODO\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* parameters (theta), matlab: second argument */\n');
    fprintf(fid, '    np = 1; // TODO\n');
    fprintf(fid, '    p = malloc(sizeof(realtype) * np);\n');
    fprintf(fid, '    // TODO p[0] = 1;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* plist, matlab: fifth argument */\n');
    fprintf(fid, '    // parameter ordering\n');
    fprintf(fid, '    plist = malloc(np * sizeof(int));\n');
    fprintf(fid, '    for (int i = 0; i < np; i++) {\n');
    fprintf(fid, '        plist[i] = i;\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* constants, matlab: third argument */\n');
    fprintf(fid, '    const int numK = 0; // TODO\n');
    fprintf(fid, '    k = malloc(sizeof(realtype) * numK);\n');
    fprintf(fid, '    // TODO k[0] = 0.1;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* Options ; matlab: fourth argument   */\n');
    fprintf(fid, '    nmaxevent = 0; // ?\n');
    fprintf(fid, '    z2event = malloc(sizeof(realtype) * ne);\n');
    fprintf(fid, '    for(int i = 0; i < ne; ++i)\n');
    fprintf(fid, '        z2event[i] = i;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    tstart = 0;\n');
    fprintf(fid, '    atol = 1E-16;\n');
    fprintf(fid, '    rtol = 1E-8;\n');
    fprintf(fid, '    maxsteps = 1e4;\n');
    fprintf(fid, '    lmm = CV_BDF;\n');
    fprintf(fid, '    iter = CV_NEWTON;\n');
    fprintf(fid, '    interpType = CV_POLYNOMIAL; //?\n');
    fprintf(fid, '    linsol = 9;\n');
    fprintf(fid, '    stldet = TRUE;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    idlist = malloc(sizeof(realtype) * np);\n');
    fprintf(fid, '    for(int i = 0; i < np; ++i)\n');
    fprintf(fid, '        idlist[i] = 0;\n');
    fprintf(fid, '    qpositivex = malloc(sizeof(realtype) * nx);\n');
    fprintf(fid, '    zeros(qpositivex, nx);\n');
    fprintf(fid, '    sensi = 0;\n');
    fprintf(fid, '    ism = 1;\n');
    fprintf(fid, '    sensi_meth = 1;\n');
    fprintf(fid, '    ordering = 0;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    //user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters] default is sensitivity initialisation based on the derivative of the state initialisation\n');
    fprintf(fid, '    b_sx0 = FALSE;\n');
    fprintf(fid, '    sx0data = 0;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* pbarm parameterscales ; matlab: sixth argument*/\n');
    fprintf(fid, '    pbar = malloc(sizeof(realtype) * np);\n');
    fprintf(fid, '    ones(pbar, np);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    //    /* xscale, matlab: seventh argument */\n');
    fprintf(fid, '    //    xbar = mxGetPr(prhs[6]);\n');
    fprintf(fid, '    xbar = 0; // xscale\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* END SET USER DATA */\n');
    fprintf(fid, '\n');
    fprintf(fid, '    processUserData(udata);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    return(udata);\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'void processUserData(UserData udata) {\n');
    fprintf(fid, '    if (nx>0) {\n');
    fprintf(fid, '        /* initialise temporary jacobian storage */\n');
    fprintf(fid, '        tmp_J = NewSparseMat(nx,nx,nnz);\n');
    fprintf(fid, '        M_tmp = malloc(nx*nx*sizeof(realtype));\n');
    fprintf(fid, '        dfdx_tmp = malloc(nx*nx*sizeof(realtype));\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '    if (sensi>0) {\n');
    fprintf(fid, '        /* initialise temporary dxdotdp storage */\n');
    fprintf(fid, '        tmp_dxdotdp = malloc(nx*np*sizeof(realtype));\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '    if (ne>0) {\n');
    fprintf(fid, '        /* initialise temporary stau storage */\n');
    fprintf(fid, '        stau_tmp = malloc(np*sizeof(realtype));\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '    w_tmp = malloc(nw*sizeof(realtype));\n');
    fprintf(fid, '    dwdx_tmp = malloc(ndwdx*sizeof(realtype));\n');
    fprintf(fid, '    dwdp_tmp = malloc(ndwdp*sizeof(realtype));\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '    udata->am_nan_dxdotdp = FALSE;\n');
    fprintf(fid, '    udata->am_nan_J = FALSE;\n');
    fprintf(fid, '    udata->am_nan_JSparse = FALSE;\n');
    fprintf(fid, '    udata->am_nan_xdot = FALSE;\n');
    fprintf(fid, '    udata->am_nan_xBdot = FALSE;\n');
    fprintf(fid, '    udata->am_nan_qBdot = FALSE;\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'ExpData getExperimentalData(UserData udata) {\n');
    fprintf(fid, '    ExpData edata = (ExpData) malloc(sizeof *edata);\n');
    fprintf(fid, '    if (edata == NULL) {\n');
    fprintf(fid, '        return(NULL);\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    // observation np * nt\n');
    fprintf(fid, '    my = malloc(sizeof(realtype) * nt * nx);\n');
    fprintf(fid, '    fillArray(my, nt * nx, NAN);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    // stdev of Y\n');
    fprintf(fid, '    ysigma =  malloc(sizeof(realtype) * nt * nx);\n');
    fprintf(fid, '    fillArray(ysigma, nt * nx, NAN);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    // event observations\n');
    fprintf(fid, '    mz = malloc(sizeof(realtype) * nz * nz);\n');
    fprintf(fid, '    fillArray(mz, nz * nz, NAN);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    zsigma =  malloc(sizeof(realtype) * nz * nz);\n');
    fprintf(fid, '    fillArray(zsigma, nz * nz, NAN);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    return(edata);\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'ReturnData getSimulationResults(UserData udata, ExpData edata) {\n');
    fprintf(fid, '    // TODO cleanup\n');
    fprintf(fid, '\n');
    fprintf(fid, '    TempData tdata; /* temporary data */\n');
    fprintf(fid, '    tdata = (TempData) malloc(sizeof *tdata);\n');
    fprintf(fid, '    if (tdata == NULL) goto freturn;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    void *ami_mem; /* pointer to cvodes memory block */\n');
    fprintf(fid, '    int status; /* general status flag */\n');
    fprintf(fid, '    if (nx>0) {\n');
    fprintf(fid, '        ami_mem = setupAMI(&status, udata, tdata);\n');
    fprintf(fid, '        if (ami_mem == NULL) goto freturn;\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    double *pstatus; /* return status flag */\n');
    fprintf(fid, '    pstatus = malloc(sizeof(double));\n');
    fprintf(fid, '\n');
    fprintf(fid, '    ReturnData rdata = initReturnData(udata, pstatus);\n');
    fprintf(fid, '    if (rdata == NULL) goto freturn;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    if (nx>0) {\n');
    fprintf(fid, '        if (edata == NULL) goto freturn;\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    int iroot = 0;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    int problem = workForwardProblem(udata, tdata, rdata, edata, &status, ami_mem, &iroot);\n');
    fprintf(fid, '    if(problem)\n');
    fprintf(fid, '        goto freturn;\n');
    fprintf(fid, '\n');
    fprintf(fid, '    booleantype setupBdone = false;\n');
    fprintf(fid, '    problem = workBackwardProblem(udata, tdata, rdata, edata, &status, ami_mem, &iroot, &setupBdone);\n');
    fprintf(fid, '    if(problem)\n');
    fprintf(fid, '        goto freturn;\n');
    fprintf(fid, '\n');
    fprintf(fid, 'freturn:\n');
    fprintf(fid, '    storeJacobianAndDerivativeInReturnData(udata, tdata, rdata);\n');
    fprintf(fid, '    freeTempDataAmiMem(udata, tdata, ami_mem, setupBdone, pstatus, status);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    return rdata;\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'ReturnData initReturnData(void *user_data, double *pstatus) {\n');
    fprintf(fid, '    //   rdata = setupReturnData(plhs, udata, pstatus);\n');
    fprintf(fid, '    /**\n');
    fprintf(fid, '     * setupReturnData initialises the return data struct\n');
    fprintf(fid, '     * @param[in] plhs user input @type mxArray\n');
    fprintf(fid, '     * @param[in] user_data pointer to the user data struct @type UserData\n');
    fprintf(fid, '     * @param[out] pstatus pointer to the flag indicating the execution status @type double\n');
    fprintf(fid, '     * @return rdata: return data struct @type ReturnData\n');
    fprintf(fid, '     */\n');
    fprintf(fid, '    ReturnData rdata; /* returned rdata struct */\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* Return rdata structure */\n');
    fprintf(fid, '    rdata = (ReturnData) malloc(sizeof *rdata);\n');
    fprintf(fid, '    if (rdata == NULL) return(NULL);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    initUserDataFields(user_data, rdata, pstatus); // TODO need to remove initfields...\n');
    fprintf(fid, '\n');
    fprintf(fid, '    return(rdata);\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'void freeTempDataAmiMem(UserData udata, TempData tdata, void *ami_mem, booleantype setupBdone, double *pstatus, int status) {\n');
    fprintf(fid, '\n');
    fprintf(fid, '    /* Free memory */\n');
    fprintf(fid, '    if(nx>0) {\n');
    fprintf(fid, '        N_VDestroy_Serial(x);\n');
    fprintf(fid, '        N_VDestroy_Serial(dx);\n');
    fprintf(fid, '        N_VDestroy_Serial(xdot);\n');
    fprintf(fid, '        N_VDestroy_Serial(x_old);\n');
    fprintf(fid, '        N_VDestroy_Serial(dx_old);\n');
    fprintf(fid, '        N_VDestroy_Serial(xdot_old);\n');
    fprintf(fid, '        DestroyMat(Jtmp);\n');
    fprintf(fid, '        DestroySparseMat(tmp_J);\n');
    fprintf(fid, '        if (ne>0) {\n');
    fprintf(fid, '            if(ami_mem) free(rootsfound);\n');
    fprintf(fid, '            if(ami_mem) free(rootvals);\n');
    fprintf(fid, '            if(ami_mem) free(rootidx);\n');
    fprintf(fid, '            if(ami_mem)    free(sigma_z);\n');
    fprintf(fid, '            if(ami_mem)    free(nroots);\n');
    fprintf(fid, '            if(ami_mem)    free(discs);\n');
    fprintf(fid, '            if(ami_mem)    free(h);\n');
    fprintf(fid, '\n');
    fprintf(fid, '            if(ami_mem)    free(deltax);\n');
    fprintf(fid, '            if(ami_mem)    free(deltasx);\n');
    fprintf(fid, '            if(ami_mem)    free(deltaxB);\n');
    fprintf(fid, '            if(ami_mem)    free(deltaqB);\n');
    fprintf(fid, '        }\n');
    fprintf(fid, '\n');
    fprintf(fid, '        if(ny>0) {\n');
    fprintf(fid, '            if(sigma_y)    free(sigma_y);\n');
    fprintf(fid, '        }\n');
    fprintf(fid, '        if (sensi >= 1) {\n');
    fprintf(fid, '            if(ami_mem)    free(dydx);\n');
    fprintf(fid, '            if(ami_mem)    free(dydp);\n');
    fprintf(fid, '            if (sensi_meth == AMI_FSA) {\n');
    fprintf(fid, '                N_VDestroyVectorArray_Serial(NVsx,np);\n');
    fprintf(fid, '            }\n');
    fprintf(fid, '            if (sensi_meth == AMI_ASA) {\n');
    fprintf(fid, '                if(status == 0) {\n');
    fprintf(fid, '                    N_VDestroyVectorArray_Serial(NVsx,np);\n');
    fprintf(fid, '                }\n');
    fprintf(fid, '            }\n');
    fprintf(fid, '\n');
    fprintf(fid, '            if (sensi_meth == AMI_FSA) {\n');
    fprintf(fid, '                N_VDestroyVectorArray_Serial(sdx, np);\n');
    fprintf(fid, '            }\n');
    fprintf(fid, '            if (sensi_meth == AMI_ASA) {\n');
    fprintf(fid, '                if(ami_mem)    free(dgdp);\n');
    fprintf(fid, '                if(ami_mem)    free(dgdx);\n');
    fprintf(fid, '                if(ami_mem)    free(drdp);\n');
    fprintf(fid, '                if(ami_mem)    free(drdx);\n');
    fprintf(fid, '                if (ne>0) {\n');
    fprintf(fid, '                    if(ami_mem)    free(dzdp);\n');
    fprintf(fid, '                    if(ami_mem)    free(dzdx);\n');
    fprintf(fid, '                }\n');
    fprintf(fid, '                if(ami_mem)     free(llhS0);\n');
    fprintf(fid, '                if(ami_mem)    free(dsigma_ydp);\n');
    fprintf(fid, '                if (ne>0) {\n');
    fprintf(fid, '                    if(ami_mem)    free(dsigma_zdp);\n');
    fprintf(fid, '                }\n');
    fprintf(fid, '                if(setupBdone)      N_VDestroy_Serial(dxB);\n');
    fprintf(fid, '                if(setupBdone)      N_VDestroy_Serial(xB);\n');
    fprintf(fid, '                if(setupBdone)      N_VDestroy_Serial(xB_old);\n');
    fprintf(fid, '                if(setupBdone)      N_VDestroy_Serial(xQB);\n');
    fprintf(fid, '                if(setupBdone)      N_VDestroy_Serial(xQB_old);\n');
    fprintf(fid, '            }\n');
    fprintf(fid, '\n');
    fprintf(fid, '        }\n');
    fprintf(fid, '        if(ami_mem)     N_VDestroy_Serial(id);\n');
    fprintf(fid, '        if(ami_mem)     AMIFree(&ami_mem);\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    if(udata)   free(plist);\n');
    fprintf(fid, '    if (sensi >= 1) {\n');
    fprintf(fid, '        if(udata)   free(tmp_dxdotdp);\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    if(tdata)    free(tdata);\n');
    fprintf(fid, '\n');
    fprintf(fid, '    if(pstatus) *pstatus = (double) status;\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'void processReturnData(ReturnData rdata, UserData udata) {\n');
    fprintf(fid, '    printReturnData(rdata, udata);\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'void printReturnData(ReturnData rdata, UserData udata) {\n');
    fprintf(fid, '\n');
    fprintf(fid, '    printf("tsdata: ");\n');
    fprintf(fid, '    printArray(tsdata, nt);\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '    printf("\\n\\nxdata\\n");\n');
    fprintf(fid, '    for(int i = 0; i < nx; ++i) {\n');
    fprintf(fid, '        for(int j = 0; j < nt; ++j)\n');
    fprintf(fid, '            printf("%%e\\t", rdata->am_xdata[j +  nt * i]);\n');
    fprintf(fid, '        printf("\\n");\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    printf("\\nydata\\n");\n');
    fprintf(fid, '    for(int i = 0; i < ny; ++i) {\n');
    fprintf(fid, '        for(int j = 0; j < nt; ++j)\n');
    fprintf(fid, '            printf("%%e\\t", rdata->am_ydata[j +  nt * i]);\n');
    fprintf(fid, '        printf("\\n");\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    printf("\\n\\nxdotdata: ");\n');
    fprintf(fid, '    for(int i = 0; i < nx; ++i)\n');
    fprintf(fid, '        printf("%%e\\t", rdata->am_xdotdata[i]);\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '    printf("\\njdata\\n");\n');
    fprintf(fid, '    for(int i = 0; i < nx; ++i) {\n');
    fprintf(fid, '        for(int j = 0; j < nx; ++j)\n');
    fprintf(fid, '            printf("%%e\\t", rdata->am_Jdata[i + nx * j]);\n');
    fprintf(fid, '        printf("\\n");\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '\n');
    fprintf(fid, '    printf("\\nnumsteps: \\t\\t");\n');
    fprintf(fid, '    printfArray(numstepsdata, nt, "%%.0f ");\n');
    fprintf(fid, '    printf("\\nnumrhsevalsdata: \\t");\n');
    fprintf(fid, '    printfArray(numrhsevalsdata, nt, "%%.0f ");\n');
    fprintf(fid, '    printf("\\norder: \\t\\t");\n');
    fprintf(fid, '    printfArray(orderdata, nt, "%%.0f ");\n');
    fprintf(fid, '\n');
    fprintf(fid, '    printf("\\n");\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'void freeExpData(ExpData edata) {\n');
    fprintf(fid, '    if(edata) {\n');
    fprintf(fid, '        free(my);\n');
    fprintf(fid, '        free(ysigma);\n');
    fprintf(fid, '        free(mz);\n');
    fprintf(fid, '        free(zsigma);\n');
    fprintf(fid, '        free(edata);\n');
    fprintf(fid, '    }\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'void freeReturnData(ReturnData rdata) {\n');
    fprintf(fid, '    // TODO ... see initReturnData\n');
    fprintf(fid, '    free(rdata);\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\n');
    fprintf(fid, 'void freeUserData(UserData udata) {\n');
    fprintf(fid, '    if(udata)\n');
    fprintf(fid, '        free(udata);\n');
    fprintf(fid, '}\n');

end
