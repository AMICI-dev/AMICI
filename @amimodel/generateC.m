function generateC(this)
% generateC generates the c files which will be used in the compilation.
%
% Return values:
%  this: model definition object @type amimodel

% different signatures for cvodes / idas
if(strcmp(this.wtype,'iw'))
    dxvec = 'dx,';
    rtcj = 'cj,';
else
    dxvec = 'NULL,';
    rtcj = '';
end


% write fun ccode

for ifun = this.funs
    if(isfield(this.fun,ifun{1}))
        fprintf([ifun{1} ' | ']);
        fid = fopen(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' ifun{1} '.cpp']),'w');
        fprintf(fid,'\n');
        fprintf(fid,'#include <include/symbolic_functions.h>\n');
        fprintf(fid,'#include <include/amici.h>\n');
        fprintf(fid,'#include <include/amici_model.h>\n');
        fprintf(fid,'#include <string.h>\n');
        if( strfind(this.fun.(ifun{1}).argstr,'user_data') )
            fprintf(fid,'#include <include/tdata.h>\n');
            fprintf(fid,'#include <include/udata.h>\n');
        end
        if( strfind(this.fun.(ifun{1}).argstr,'tdata') )
            fprintf(fid,'#include <include/tdata.h>\n');
            fprintf(fid,'#include <include/udata.h>\n');
        end
        if( strfind(this.fun.(ifun{1}).argstr,'rdata') )
            fprintf(fid,'#include <include/rdata.h>\n');
        end
        if( strfind(this.fun.(ifun{1}).argstr,'edata') )
            fprintf(fid,'#include <include/edata.h>\n');
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
        elseif( strcmp(ifun{1},'J') || strcmp(ifun{1},'JDiag') || strcmp(ifun{1},'JB') || strcmp(ifun{1},'JSparse') || strcmp(ifun{1},'JSparseB') || strcmp(ifun{1},'xBdot') || strcmp(ifun{1},'dfdx'))
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
                fprintf(fid,'TempData *tdata = (TempData*) user_data;\n');
                fprintf(fid,'Model *model = (Model*) tdata->model;\n');
                fprintf(fid,'UserData *udata = (UserData*) tdata->udata;\n');
            end
            if( strfind(this.fun.(ifun{1}).argstr,'tdata') )
                fprintf(fid,'Model *model = (Model*) tdata->model;\n');
                fprintf(fid,'UserData *udata = (UserData*) tdata->udata;\n');
            end
            this.fun.(ifun{1}).printLocalVars(this,fid);
            if(~isempty(strfind(this.fun.(ifun{1}).argstr,'N_Vector x')) && ~isempty(strfind(this.fun.(ifun{1}).argstr,'realtype t')))
                if(or(not(strcmp(this.wtype,'iw')),~isempty(strfind(this.fun.(ifun{1}).argstr,'N_Vector dx'))))
                    if(~any(ismember(ifun{1},{'w','sxdot','dxdotdp','dfdx','qBdot'})))
                        fprintf(fid,['status = w_' this.modelname '(t,x,' dxvec 'tdata);\n']);
                    end
                end
            end
            if( strcmp(ifun{1},'sxdot') )
                if(strcmp(this.wtype,'iw'))
                    fprintf(fid,'int ip;\n');
                    fprintf(fid,['status = dfdx_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                    fprintf(fid,['status = M_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                    fprintf(fid,['status = dxdotdp_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                    fprintf(fid,'for(ip = 0; ip<udata->nplist; ip++) {\n');
                    fprintf(fid,'sx_tmp = N_VGetArrayPointer(sx[udata->plist[ip]]);\n');
                    fprintf(fid,'sdx_tmp = N_VGetArrayPointer(sdx[udata->plist[ip]]);\n');
                    fprintf(fid,'sxdot_tmp = N_VGetArrayPointer(sxdot[udata->plist[ip]]);\n');
                    fprintf(fid,['memset(sxdot_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                    this.fun.(ifun{1}).writeCcode(this,fid);
                    fprintf(fid,'}\n');
                else
                    fprintf(fid,'if(ip == 0) {\n');
                    fprintf(fid,['    status = JSparse_' this.modelname '(t,' rtcj 'x,xdot,tdata->J,user_data,NULL,NULL,NULL);\n']);
                    fprintf(fid,['    status = dxdotdp_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                    fprintf(fid,'}\n');
                    this.fun.(ifun{1}).writeCcode(this,fid);
                end
            elseif( strcmp(ifun{1},'qBdot') )
                fprintf(fid,['status = dwdp_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                fprintf(fid,'for(ip = 0; ip<udata->nplist; ip++) {\n');
                fprintf(fid,'switch (udata->plist[ip]) {\n');
                this.fun.(ifun{1}).writeCcode_sensi(this,fid);
                fprintf(fid,'}\n');
                fprintf(fid,'}\n');
            elseif(this.fun.(ifun{1}).sensiflag)
                if( strcmp(ifun{1},'dxdotdp'))
                    fprintf(fid,['status = dwdp_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                end
                fprintf(fid,'for(ip = 0; ip<udata->nplist; ip++) {\n');
                if(ismember('*sx',this.fun.(ifun{1}).nvecs))
                    fprintf(fid,'sx_tmp = N_VGetArrayPointer(sx[ip]);\n');
                end
                if(ismember('*sx0',this.fun.(ifun{1}).nvecs))
                    fprintf(fid,'sx0_tmp = N_VGetArrayPointer(sx0[ip]);\n');
                    fprintf(fid,['memset(sx0_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                end
                if(ismember('*sdx0',this.fun.(ifun{1}).nvecs))
                    fprintf(fid,'sdx0_tmp = N_VGetArrayPointer(sdx0[ip]);\n');
                    fprintf(fid,['memset(sdx0_tmp,0,sizeof(realtype)*' num2str(this.nx) ');\n']);
                end
                fprintf(fid,'switch (udata->plist[ip]) {\n');
                this.fun.(ifun{1}).writeCcode_sensi(this,fid);
                fprintf(fid,'}\n');
                fprintf(fid,'}\n');
            else
                if( strcmp(ifun{1},'J') || strcmp(ifun{1},'JDiag') || strcmp(ifun{1},'JB') || strcmp(ifun{1},'JSparse') || strcmp(ifun{1},'JSparseB') || strcmp(ifun{1},'xBdot') || strcmp(ifun{1},'dfdx') )
                    fprintf(fid,['status = dwdx_' this.modelname '(t,x,' dxvec 'user_data);\n']);
                end
                this.fun.(ifun{1}).writeCcode(this,fid);
            end
            if(strcmp(ifun{1},'dxdotdp'))
                fprintf(fid,'for(ip = 0; ip<udata->nplist; ip++) {\n');
                fprintf(fid,'   for(ix = 0; ix<model->nx; ix++) {\n');
                fprintf(fid,'       if(amiIsNaN(tdata->dxdotdp[ix+ip*model->nx])) {\n');
                fprintf(fid,'           tdata->dxdotdp[ix+ip*model->nx] = 0;\n');
                fprintf(fid,'           if(!tdata->nan_dxdotdp) {\n');
                fprintf(fid,'               warnMsgIdAndTxt("AMICI:mex:fdxdotdp:NaN","AMICI replaced a NaN value in dxdotdp and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'               tdata->nan_dxdotdp = TRUE;\n');
                fprintf(fid,'           }\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'       if(amiIsInf(tdata->dxdotdp[ix+ip*model->nx])) {\n');
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
                    fprintf(fid,'       if(!tdata->nan_xdot) {\n');
                    fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fxdot:NaN","AMICI replaced a NaN value in xdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                    fprintf(fid,'           tdata->nan_xdot = TRUE;\n');
                    fprintf(fid,'       }\n');
                    fprintf(fid,'   }\n');
                    fprintf(fid,'   if(amiIsInf(xdot_tmp[ix])) {\n');
                    fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fxdot:Inf","AMICI encountered an Inf value in xdot! Aborting simulation ... ");\n');
                    fprintf(fid,'       return(-1);\n');
                    fprintf(fid,'   }');
                    fprintf(fid,'   if(udata->qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) {\n');
                    fprintf(fid,'       xdot_tmp[ix] = -xdot_tmp[ix];\n');
                    fprintf(fid,'   }\n');
                    fprintf(fid,'}\n');
                end
            end
            if(strcmp(ifun{1},'JSparse'))
                fprintf(fid,['for(inz = 0; inz<' num2str(this.nnz) '; inz++) {\n']);
                fprintf(fid,'   if(amiIsNaN(J->data[inz])) {\n');
                fprintf(fid,'       J->data[inz] = 0;\n');
                fprintf(fid,'       if(!tdata->nan_JSparse) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           tdata->nan_JSparse = TRUE;\n');
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
                fprintf(fid,'       if(!tdata->nan_J) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fJ:NaN","AMICI replaced a NaN value in Jacobian and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           tdata->nan_J = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'   if(amiIsInf(J->data[ix])) {\n');
                fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fJ:Inf","AMICI encountered an Inf value in Jacobian! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'JDiag'))
                fprintf(fid,['for(ix = 0; ix<' num2str(this.nx) '; ix++) {\n']);
                fprintf(fid,'   if(amiIsNaN(JDiag_tmp[ix])) {\n');
                fprintf(fid,'       JDiag_tmp[ix] = 0;\n');
                fprintf(fid,'       if(!tdata->nan_JDiag) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fJDiag:NaN","AMICI replaced a NaN value on Jacobian diagonal and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           tdata->nan_JDiag = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'   if(amiIsInf(JDiag_tmp[ix])) {\n');
                fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fJDiag:Inf","AMICI encountered an Inf value on Jacobian diagonal! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }\n');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'xBdot'))
                fprintf(fid,['for(ix = 0; ix<' num2str(this.nx) '; ix++) {\n']);
                fprintf(fid,'   if(amiIsNaN(xBdot_tmp[ix])) {\n');
                fprintf(fid,'       xBdot_tmp[ix] = 0;');
                fprintf(fid,'       if(!tdata->nan_xBdot) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fxBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           tdata->nan_xBdot = TRUE;\n');
                fprintf(fid,'       }\n');
                fprintf(fid,'   }');
                fprintf(fid,'   if(amiIsInf(xBdot_tmp[ix])) {\n');
                fprintf(fid,'       warnMsgIdAndTxt("AMICI:mex:fxBdot:Inf","AMICI encountered an Inf value in xBdot! Aborting simulation ... ");\n');
                fprintf(fid,'       return(-1);\n');
                fprintf(fid,'   }');
                fprintf(fid,'}\n');
            end
            if(strcmp(ifun{1},'qBdot'))
                fprintf(fid,'for(ip = 0; ip<udata->nplist*model->nJ; ip++) {\n');
                fprintf(fid,'   if(amiIsNaN(qBdot_tmp[ip])) {\n');
                fprintf(fid,'       qBdot_tmp[ip] = 0;');
                fprintf(fid,'       if(!tdata->nan_qBdot) {\n');
                fprintf(fid,'           warnMsgIdAndTxt("AMICI:mex:fqBdot:NaN","AMICI replaced a NaN value in xBdot and replaced it by 0.0. This will not be reported again for this simulation run.");\n');
                fprintf(fid,'           tdata->nan_qBdot = TRUE;\n');
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
% function declarations for this.modelname.cpp
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
fprintf(fid,['#endif /* _am_' this.modelname '_h */\n']);
fclose(fid);

%
%----------------------------------------------------------------
% funname.h
% function declarations for funname.cpp
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
    fprintf(fid,'#include <sundials/sundials_types.h>\n');
    fprintf(fid,'#include <sundials/sundials_nvector.h>\n');
    fprintf(fid,'#include <sundials/sundials_sparse.h>\n');
    fprintf(fid,'#include <sundials/sundials_direct.h>\n\n');
    fprintf(fid,'class UserData;\nclass ReturnData;\nclass TempData;\nclass ExpData;\n\n');
    fprintf(fid,['int ' ifun{1} '_' this.modelname '' fun.argstr ';\n']);
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['#endif /* _am_' this.modelname '_' ifun{1} '_h */\n']);
    fclose(fid);
end
%
%----------------------------------------------------------------
% wrapfunctions.cpp
% function definitions
%----------------------------------------------------------------
%

fprintf('wrapfunctions | ');
fid = fopen(fullfile(this.wrap_path,'models',this.modelname,'wrapfunctions.cpp'),'w');

fprintf(fid,'#include "wrapfunctions.h"\n');
fprintf(fid,'#include <include/amici_model.h>\n');
fprintf(fid,'#include <include/udata.h>\n');

if(~strcmp(this.wtype,'iw'))
    fprintf(fid,'#include <include/cvodewrap.h>\n');
else
    fprintf(fid,'#include <include/idawrap.h>\n');
end
fprintf(fid,'\n');

fprintf(fid,'Solver *getSolver(){\n');
if(strcmp(this.wtype,'iw'))
    fprintf(fid, '    return new IDASolver();\n');
else
    fprintf(fid, '    return new CVodeSolver();\n');
end
fprintf(fid,'}\n\n');

fprintf(fid,'Model *getModel() {\n');
    fprintf(fid, ['    return new Model_' this.modelname '();\n']);
fprintf(fid,'}\n\n');

ffuns = {'x0','dx0','sx0','sdx0','J','JB','JDiag','Jv','root','rz','srz','stau',...
    'y','dydp','dydx','z','sz','dzdp','dzdx','drzdp','drzdx', 'sxdot', ...
    'xdot','xBdot','qBdot','dxdotdp','deltax','deltasx','deltaxB','deltaqB',...
    'sigma_y','dsigma_ydp','sigma_z','dsigma_zdp',...
    'JSparse', 'JBand', 'JSparseB', 'JBandB', 'JvB', ...
    'Jy','Jz','Jrz','dJydy','dJydsigma','dJzdz','dJzdsigma','dJrzdz','dJrzdsigma'};

for iffun = ffuns
    % check whether the function was generated, otherwise generate (but
    % whithout symbolic expressions)
    if(isfield(this.fun,iffun{1}))
        fun = this.fun.(iffun{1});
    else
        fun = amifun(iffun{1},this);
    end
    fprintf(fid,['int f' iffun{1} fun.fargstr '{\n']);
    % if the function was generated, we can return it, otherwise return
    % an error
    if(ismember(iffun{1},this.funs))
        fprintf(fid,['    return ' iffun{1} '_' this.modelname removeTypes(fun.argstr) ';\n']);
    else
        if(strcmp(iffun{1},'sx0') || strcmp(iffun{1},'dx0') || strcmp(iffun{1},'sdx0'))
            % these two should always work, if they are not required
            % they should act as placeholders
            fprintf(fid,'    return 0;\n');
        else
            fprintf(fid,['    warnMsgIdAndTxt("AMICI:mex:' iffun{1} ':NotAvailable","ERROR: The function ' iffun{1} ' was called but not compiled for this model.");\n']);
            fprintf(fid,'    return -1;\n');
        end
    end
    fprintf(fid,'}\n\n');
end

fclose(fid);

%
%----------------------------------------------------------------
% wrapfunctions.h
% function declarations for wrapfunctions.cpp
%----------------------------------------------------------------
%

fid = fopen(fullfile(this.wrap_path,'models',this.modelname,'wrapfunctions.h'),'w');
fprintf(fid,'#ifndef _amici_wrapfunctions_h\n');
fprintf(fid,'#define _amici_wrapfunctions_h\n');
fprintf(fid,'#include <math.h>\n');
fprintf(fid,'#include <include/amici_model.h>\n');
fprintf(fid,'\n');
fprintf(fid,['#include "' this.modelname '.h"\n']);
fprintf(fid,'\n');
fprintf(fid,'class UserData;\nclass Solver;\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'#define pi M_PI\n');
fprintf(fid,'\n');
fprintf(fid,'#ifdef __cplusplus\n#define EXTERNC extern "C"\n#else\n#define EXTERNC\n#endif\n');
fprintf(fid,'\n');
fprintf(fid,'UserData getUserData();\n');
fprintf(fid,'Solver *getSolver();\n');
fprintf(fid,'Model *getModel();\n');

for iffun = ffuns
    % check whether the function was generated, otherwise generate (but
    % whithout symbolic expressions)
    if(isfield(this.fun,iffun{1}))
        fun = this.fun.(iffun{1});
    else
        fun = amifun(iffun{1},this);
    end
    fprintf(fid,['int f' iffun{1} fun.fargstr ';\n']);
end

% Subclass Model
fprintf(fid,'\n');
fprintf(fid,['class Model_' this.modelname ' : public Model {\n']);
fprintf(fid,'public:\n');
fprintf(fid,['    Model_' this.modelname '() : Model(' num2str(this.np) ',\n']);
fprintf(fid,['                    ' num2str(this.nx) ',\n']);
fprintf(fid,['                    ' num2str(this.nxtrue) ',\n']);
fprintf(fid,['                    ' num2str(this.nk) ',\n']);
fprintf(fid,['                    ' num2str(this.ny) ',\n']);
fprintf(fid,['                    ' num2str(this.nytrue) ',\n']);
fprintf(fid,['                    ' num2str(this.nz) ',\n']);
fprintf(fid,['                    ' num2str(this.nztrue) ',\n']);
fprintf(fid,['                    ' num2str(this.nevent) ',\n']);
fprintf(fid,['                    ' num2str(this.ng) ',\n']);
fprintf(fid,['                    ' num2str(this.nw) ',\n']);
fprintf(fid,['                    ' num2str(this.ndwdx) ',\n']);
fprintf(fid,['                    ' num2str(this.ndwdp) ',\n']);
fprintf(fid,['                    ' num2str(this.nnz) ',\n']);
fprintf(fid,['                    ' num2str(this.ubw) ',\n']);
fprintf(fid,['                    ' num2str(this.lbw) ',\n']);
switch(this.o2flag)
    case 1
        fprintf(fid,'                    AMICI_O2MODE_FULL)\n');
    case 2
        fprintf(fid,'                    AMICI_O2MODE_DIR)\n');
    otherwise
        fprintf(fid,'                    AMICI_O2MODE_NONE)\n');
end
fprintf(fid,'{\n');
fprintf(fid,['    z2event = new int[nz] {' num2str(transpose(this.z2event), '%d, ') '};\n']);
fprintf(fid,['    idlist = new realtype[nx] {' num2str(transpose(double(this.id)), '%d, ') '};\n']);
fprintf(fid,'}\n\n');

for iffun = this.funs
    % check whether the function was generated, otherwise generate (but
    % whithout symbolic expressions)
    if(isfield(this.fun,iffun{1}))
        fun = this.fun.(iffun{1});
    else
        fun = amifun(iffun{1},this);
    end
    fprintf(fid,['    int f' iffun{1} fun.fargstr ' {\n']);
    fprintf(fid,['        return ' iffun{1} '_' this.modelname removeTypes(fun.argstr) ';\n']);
    fprintf(fid,'    }\n\n');
end
fprintf(fid,'};\n\n');

fprintf(fid,'#endif /* _amici_wrapfunctions_h */\n');
fclose(fid);

fprintf('CMakeLists | ');
generateCMakeFile(this);

fprintf('main | ');
generateMainC(this);

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
argstr = strrep(argstr,'TempData','');
argstr = strrep(argstr,'ReturnData','');
argstr = strrep(argstr,'const ExpData','');
argstr = strrep(argstr,'N_Vector','');
argstr = strrep(argstr,'long','');
argstr = strrep(argstr,'DlsMat','');
argstr = strrep(argstr,'SlsMat','');
argstr = strrep(argstr,'*','');
argstr = strrep(argstr,' ','');
argstr = strrep(argstr,',',', ');

end


function generateCMakeFile(this)
    sourceStr = '';
    for j=1:length(this.funs)
        funcName = this.funs{j};
        sourceStr = [ sourceStr, sprintf('${MODEL_DIR}/%s_%s.cpp\n', this.modelname, funcName) ];
    end
    
    t = template();
    t.add('TPL_MODELNAME', this.modelname);
    t.add('TPL_SOURCES', sourceStr);
    CMakeFileName = fullfile(this.wrap_path,'models',this.modelname,'CMakeLists.txt');
    CMakeTemplateFileName = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src/CMakeLists.template.txt');
    t.replace(CMakeTemplateFileName, CMakeFileName);
end

function generateMainC(this)
    mainFileSource = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src/main.template.cpp');
    mainFileDestination = fullfile(this.wrap_path,'models',this.modelname,'main.cpp');
    copyfile(mainFileSource, mainFileDestination);
end
