% function [] = generateAMICIfile(model,name,path,atol,rtol,maxsteps,compile_forward,compile_adjoint)
function [] = generateAMICIfile(model,name,path,compile_forward,compile_adjoint)

% Open file
fid = fopen([path name '_syms.m'],'w');

% Function header
% str = ['function [model] = ' name '_syms() \n\n'...
%        '%%%% CVODES OPTIONS\n'...
%        '%% set the default absolute tolerance\n'...
%        'model.atol = ' num2str(atol,'%10.5e') ';\n'...
%        '%% set the default relative tolerance\n'...
%        'model.rtol = ' num2str(rtol,'%10.5e') ';\n'...
%        '%% set the default maximum number of integration steps\n'...
%        'model.maxsteps = ' num2str(maxsteps,'%i') ';\n'...
str = ['function [model] = ' name '_syms() \n\n'...
       '%% set the parametrisation of the problem options are ''log'', ''log10'' and ''lin'' (default)\n'...
       'model.param = ''log10'';\n'...
       '%% set compiler options\n'...
       'model.forward = ' compile_forward ';\n'...
       'model.adjoint = ' compile_adjoint ';\n\n'];

% State variables
str = [str '%%%% STATE VARIABLES\n'...
           '%% create state variables syms\n'];
for i = 1:length(model.sym.x)
    str = [str 'syms ' char(model.sym.x(i)) '\n'];
end
str = [str '\n' '%% create state variable vector\n'...
                'x = ['];
for i = 1:length(model.sym.x)
    if i >= 2
        str = [str '     '];
    end
    str = [str char(model.sym.x(i))];
    if i < length(model.sym.x)
        str = [str '\n'];
    end
end
str = [str '];\n\n'];

% Parameters
str = [str '%%%% PARAMETERS\n'...
           '%% create parameters syms\n'];
for i = 1:length(model.sym.p)
    str = [str 'syms ' char(model.sym.p(i)) '\n'];
end
str = [str '\n' '%% create parameter vector\n'...
                'p = ['];
for i = 1:length(model.sym.p)
    if i >= 2
        str = [str '     '];
    end
    str = [str char(model.sym.p(i))];
    if i < length(model.sym.p)
        str = [str '\n'];
    end
end
str = [str '];\n\n'];

% Inputs
str = [str '%%%% CONSTANTS\n'...
           '%% create constants syms\n'];
for i = 1:length(model.sym.k)
    str = [str 'syms ' char(model.sym.k(i)) '\n'];
end
str = [str '\n' '%% create constant vector\n'...
                'k = ['];
for i = 1:length(model.sym.k)
    if i >= 2
        str = [str '     '];
    end
    str = [str char(model.sym.k(i))];
    if i < length(model.sym.k)
        str = [str '\n'];
    end
end
str = [str '];\n\n'];

% Right-hand side of differential equation
str = [str '%%%% RIGHT-HAND SIDE OF DIFFERENTIAL EQUATION\n'...
           '%% create symbolic variable for time\n'...
           'syms t\n\n'...
           '%% symbolic expression for right-hand side of differential equation\n'...
           'xdot = sym(zeros(' num2str(length(model.sym.x)) ',1));\n'];
for i = 1:length(model.sym.xdot)
    str = [str 'xdot(' num2str(i,'%i') ') = ' char(model.sym.xdot(i)) ';\n'];
end
str = [str '\n'];

% Initial condition
str = [str '%%%% INITIAL CONDITION\n'...
           'x0 = sym(zeros(' num2str(length(model.sym.x0),'%i') ',1));\n\n'];
for i = 1:length(model.sym.x0)
    if isnumeric(model.sym.x0(i))
        str = [str 'x0(' num2str(i,'%i') ') = ' num2str(model.sym.x0(i)) ';\n'];
    else
        str = [str 'x0(' num2str(i,'%i') ') = ' char(model.sym.x0(i)) ';\n'];
    end
end
str = [str '\n'];

% Observable
str = [str '%%%% OBSERVABLE\n'...
           'y = sym(zeros(' num2str(length(model.sym.y),'%i') ',1));\n\n'];
for i = 1:length(model.sym.y)
    str = [str 'y(' num2str(i,'%i') ') = ' char(model.sym.y(i)) ';\n'];
end
str = [str '\n'];

% Events
str = [str '%%%% EVENTS\n'...
           '%% this part is optional and can be ommited\n\n'...
           '%% events fire when there is a zero crossing of the root function\n'...
           'event = sym(zeros(' num2str(length(model.sym.event),'%i') ',1));\n\n'];

for i = 1:length(model.sym.event)
    str = [str 'event(' num2str(i,'%i') ') = ' char(model.sym.event(i)) ';\n'];
end

% System struct
str = [str '%%%% SYSTEM STRUCT\n'...
           'model.sym.x = x;\n'...
           'model.sym.k = k;\n'...
           'model.sym.event = event;\n'...
           'model.sym.xdot = xdot;\n'...
           'model.sym.p = p;\n'...
           'model.sym.x0 = x0;\n'...
           'model.sym.y = y;\n\n'...
           'end'];

% Write to file
fprintf(fid,str);
fclose(fid);