function generateC(this)
% generateC generates the c files which will be used in the compilation.
%
% Return values:
%  void

% different signatures for cvodes / idas
if(strcmp(this.wtype,'iw'))
    dxvec = 'dx,';
    rtcj = 'cj,';
else
    dxvec = 'NULL,';
    rtcj = '0.0,';
end


% write fun ccode

for ifun = this.funs
    if(isfield(this.fun,ifun{1}))
        bodyNotEmpty = any(this.fun.(ifun{1}).sym(:)~=0);
        if(strcmp(ifun{1},'JSparse'))
            bodyNotEmpty = any(this.fun.J.sym(:)~=0);
        end
        if(strcmp(ifun{1},'JSparseB'))
            bodyNotEmpty = any(this.fun.JB.sym(:)~=0);
        end
        
        if(bodyNotEmpty)
            fprintf([ifun{1} ' | ']);
            fid = fopen(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' ifun{1} '.cpp']),'w');
            fprintf(fid,'\n');
            fprintf(fid,'#include <include/symbolic_functions.h>\n');
            fprintf(fid,'#include <include/amici_defines.h> //realtype definition\n');
            
            if(ismember(ifun{1},{'JSparse','JSparseB'}))
                fprintf(fid,'#include <sundials/sundials_sparse.h> //SlsMat definition\n');
            end
            
            fprintf(fid,'typedef amici::realtype realtype;\n');
            fprintf(fid,'#include <cmath> \n');
            fprintf(fid,'\n');
            
            % function definition
            fprintf(fid,['void ' ifun{1} '_' this.modelname '' this.fun.(ifun{1}).argstr ' {\n']);
            if(strcmp(ifun{1},'JSparse'))
                for i = 1:length(this.rowvals)
                    fprintf(fid,['  JSparse->indexvals[' num2str(i-1) '] = ' num2str(this.rowvals(i)) ';\n']);
                end
                for i = 1:length(this.colptrs)
                    fprintf(fid,['  JSparse->indexptrs[' num2str(i-1) '] = ' num2str(this.colptrs(i)) ';\n']);
                end
            end
            if(strcmp(ifun{1},'JSparseB'))
                for i = 1:length(this.rowvalsB)
                    fprintf(fid,['  JSparseB->indexvals[' num2str(i-1) '] = ' num2str(this.rowvalsB(i)) ';\n']);
                end
                for i = 1:length(this.colptrsB)
                    fprintf(fid,['  JSparseB->indexptrs[' num2str(i-1) '] = ' num2str(this.colptrsB(i)) ';\n']);
                end
            end
            
            if(strcmp(ifun{1},'JBand'))
                fprintf(fid,['return(J_' this.modelname removeTypes(this.fun.J.argstr) ');']);
            elseif(strcmp(ifun{1},'JBandB'))
                fprintf(fid,['return(JB_' this.modelname removeTypes(this.fun.JB.argstr) ');']);
            else
                if( strcmp(ifun{1},'qBdot') )
                    fprintf(fid,'switch (ip) {\n');
                    this.fun.(ifun{1}).writeCcode_sensi(this,fid);
                    fprintf(fid,'}\n');
                elseif(this.fun.(ifun{1}).sensiflag)
                    fprintf(fid,'switch (ip) {\n');
                    this.fun.(ifun{1}).writeCcode_sensi(this,fid);
                    fprintf(fid,'}\n');
                else
                    this.fun.(ifun{1}).writeCcode(this,fid);
                end
            end
            fprintf(fid,'}\n');
            fprintf(fid,'\n');
            fclose(fid);
        end
    end
end

%
%----------------------------------------------------------------
% wrapfunctions.h
% function declarations for wrapfunctions.cpp
%----------------------------------------------------------------
%

fid = fopen(fullfile(this.wrap_path,'models',this.modelname,'wrapfunctions.h'),'w');
fprintf(fid,'#ifndef _amici_wrapfunctions_h\n');
fprintf(fid,'#define _amici_wrapfunctions_h\n');
fprintf(fid,'#include <cmath>\n');
fprintf(fid,'#include <memory>\n');
fprintf(fid,'#include <include/amici_defines.h>\n');
fprintf(fid,'#include <sundials/sundials_sparse.h> //SlsMat definition\n');
if(~strcmp(this.wtype,'iw'))
    fprintf(fid,'#include <include/amici_solver_cvodes.h>\n');
    fprintf(fid,'#include <include/amici_model_ode.h>\n');
else
    fprintf(fid,'#include <include/amici_solver_idas.h>\n');
    fprintf(fid,'#include <include/amici_model_dae.h>\n');
end
fprintf(fid,'\n');
fprintf(fid,'namespace amici {\nclass Solver;\n}\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'#define pi M_PI\n');
fprintf(fid,'\n');
fprintf(fid,'#ifdef __cplusplus\n#define EXTERNC extern "C"\n#else\n#define EXTERNC\n#endif\n');
fprintf(fid,'\n');
fprintf(fid,'std::unique_ptr<amici::Model> getModel();\n');

for ifun = this.funs
    if(~isfield(this.fun,ifun{1}))
        
        this.fun(1).(ifun{1}) = amifun(ifun{1},this); % don't use getfun here
        % as we do not want symbolics to be generated, we only want to be able
        % access argstr
    end
    if(checkIfFunctionBodyIsNonEmpty(this,ifun{1}))
        fprintf(fid,['extern void ' ifun{1} '_' this.modelname this.fun.(ifun{1}).argstr ';\n']);
    end
end

% Subclass Model
fprintf(fid,'\n');
if(strcmp(this.wtype,'iw'))
    baseclass = 'Model_DAE';
else
    baseclass = 'Model_ODE';
end

fprintf(fid,['class Model_' this.modelname ' : public amici::' baseclass ' {\n']);
fprintf(fid,'public:\n');
fprintf(fid,['    Model_' this.modelname '() : amici::' baseclass '(' num2str(this.nx) ',\n']);
fprintf(fid,['                    ' num2str(this.nxtrue) ',\n']);
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
        fprintf(fid,'                    amici::AMICI_O2MODE_FULL,\n');
    case 2
        fprintf(fid,'                    amici::AMICI_O2MODE_DIR,\n');
    otherwise
        fprintf(fid,'                    amici::AMICI_O2MODE_NONE,\n');
end
fprintf(fid,['                    std::vector<realtype>(' num2str(this.np) '),\n']);
fprintf(fid,['                    std::vector<realtype>(' num2str(this.nk) '),\n']);
fprintf(fid,'                    std::vector<int>(),\n');
initstr = num2str(transpose(double(this.id)), '%d, ');
fprintf(fid,['                    std::vector<realtype>{' initstr(1:end-1) '},\n']);
initstr = num2str(transpose(this.z2event), '%d, ');
fprintf(fid,['                    std::vector<int>{' initstr(1:end-1) '})\n']);
fprintf(fid,['                    {};\n\n']);

for ifun = this.funs
    fprintf(fid,['    virtual void f' ifun{1} this.fun.(ifun{1}).argstr ' override {\n']);
    if(checkIfFunctionBodyIsNonEmpty(this,ifun{1}))
        fprintf(fid,['        ' ifun{1} '_' this.modelname '' removeTypes(this.fun.(ifun{1}).argstr) ';\n']);
    end
    fprintf(fid,'    }\n\n');
end
fprintf(fid,'};\n\n');

fprintf(fid,'#endif /* _amici_wrapfunctions_h */\n');
fclose(fid);


fid = fopen(fullfile(this.wrap_path,'models',this.modelname,'wrapfunctions.cpp'),'w');
fprintf(fid,'#include <include/amici_model.h>\n');
fprintf(fid,'#include "wrapfunctions.h"\n\n');
fprintf(fid,'std::unique_ptr<amici::Model> getModel() {\n');
fprintf(fid, ['    return std::unique_ptr<amici::Model>(new Model_' this.modelname '());\n']);
fprintf(fid,'}\n\n');
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
argstr = strrep(argstr,'const','');
argstr = strrep(argstr,'double','');
argstr = strrep(argstr,'SlsMat','');
argstr = strrep(argstr,'*','');
argstr = strrep(argstr,' ','');
argstr = strrep(argstr,',',', ');

end


function generateCMakeFile(this)
    sourceStr = '';
    for j=1:length(this.funs)
        funcName = this.funs{j};
        if(checkIfFunctionBodyIsNonEmpty(this,funcName))
            sourceStr = [ sourceStr, sprintf('${MODEL_DIR}/%s_%s.cpp\n', this.modelname, funcName) ];
        end
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

function nonempty = checkIfFunctionBodyIsNonEmpty(this,ifun)
    % if we don't have symbolic variables, it might have been generated before and symbolic expressions were simply not
    % regenerated. any() for empty (no generated) variables is always false.
    nonempty = or(exist(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' ifun '.cpp']),'file'),any(this.fun.(ifun).sym(:)~=0));
end
