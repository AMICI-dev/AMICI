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
        fprintf([ifun{1} ' | ']);
        fid = fopen(fullfile(this.wrap_path,'models',this.modelname,[this.modelname '_' ifun{1} '.cpp']),'w');
        fprintf(fid,'\n');
        fprintf(fid,'#include <include/symbolic_functions.h>\n');
        fprintf(fid,['#include "' this.modelname '_w.h"\n']);
        fprintf(fid,'\n');
        
        fprintf(fid,['using namespace ' this.modelname ';\n\n']);
        
        % function definition
        fprintf(fid,['void ' ifun{1} '_' this.modelname '' this.fun.(ifun{1}).argstr ' {\n']);
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
fprintf(fid,['#include "' this.modelname '.h"\n']);
if(~strcmp(this.wtype,'iw'))
    fprintf(fid,'#include <include/amici_solver_cvodes.h>\n');
else
    fprintf(fid,'#include <include/amici_solver_idas.h>\n');
end
fprintf(fid,'\n');
fprintf(fid,'namespace amici {\nclass UserData;\nclass Solver;\n}\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'#define pi M_PI\n');
fprintf(fid,'\n');
fprintf(fid,'#ifdef __cplusplus\n#define EXTERNC extern "C"\n#else\n#define EXTERNC\n#endif\n');
fprintf(fid,'\n');
fprintf(fid,'amici::UserData getUserData();\n');
fprintf(fid,'amici::Solver *getSolver();\n');
fprintf(fid,'amici::Model *getModel();\n');

% Subclass Model
fprintf(fid,'\n');
if(strcmp(this.wtype,'iw'))
    baseclass = 'Model_DAE';
else
    baseclass = 'Model_ODE';
end


fprintf(fid,['class Model_' this.modelname ' : public amici::' baseclass ' {\n']);
fprintf(fid,'public:\n');
fprintf(fid,['    Model_' this.modelname '() : amici::' baseclass '(' num2str(this.np) ',\n']);
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
        fprintf(fid,'                    amici::AMICI_O2MODE_FULL)\n');
    case 2
        fprintf(fid,'                    amici::AMICI_O2MODE_DIR)\n');
    otherwise
        fprintf(fid,'                    amici::AMICI_O2MODE_NONE)\n');
end
fprintf(fid,'    {\n');
fprintf(fid,['        z2event = new int[' num2str(this.nz) '] {' num2str(transpose(this.z2event), '%d, ') '};\n']);
fprintf(fid,['        idlist = new realtype[' num2str(this.nx) '] {' num2str(transpose(double(this.id)), '%d, ') '};\n']);
fprintf(fid,'    }\n\n');

for ifun = this.funs
    fprintf(fid,['    void ' ifun{1} '_' this.modelname  this.fun.(ifun{1}).argstr ';\n']);
    fprintf(fid,['    model_' ifun{1} ' = &' ifun{1} '_' this.modelname ';\n\n']);
end



fprintf(fid,'    amici::Solver *getSolver() override {\n');
if(strcmp(this.wtype,'iw'))
    fprintf(fid, '        return new amici::IDASolver();\n');
else
    fprintf(fid, '        return new amici::CVodeSolver();\n');
end
fprintf(fid,'    }\n\n');

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
argstr = strrep(argstr,'amici::TempData','');
argstr = strrep(argstr,'amici::ReturnData','');
argstr = strrep(argstr,'const amici::ExpData','');
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
