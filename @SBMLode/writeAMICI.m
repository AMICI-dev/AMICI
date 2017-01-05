function writeAMICI(this,modelname)
    % writeAMICI writes the symbolic information from an SBMLode object
    % into an AMICI model definition file
    %
    % Parameters:
    %  modelname: target name of the model (_syms.m will be appended to the name )
    %
    % Return values:
    fprintf('writing file ...\n')
    fid = fopen([modelname '_syms.m'],'w');
    
    fprintf(fid,['function model = ' modelname '_syms()\n']);
    fprintf(fid,'\n');
    if(strcmp(this.time_symbol,''))
        fprintf(fid,'t = sym(''t'');\n');
    else
        fprintf(fid,[this.time_symbol ' = sym(''t'');\n']);
    end
%     fprintf(fid,'model.debug = true;\n');
    writeDefinition('STATES','x','state',this,fid)
    writeDefinition('PARAMETERS','p','parameter',this,fid)
    writeDefinition('CONDITIONS','k','condition',this,fid)
    writeDerived('DYNAMICS','xdot','xdot',this,fid)
    writeDerived('INITIALIZATION','x0','initState',this,fid)
    writeDerived('OBSERVABLES','y','observable',this,fid)
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['% EVENTS\n']);
    for ievent = 1:length(this.trigger)
        str_trigger = char(this.trigger(ievent));
        str_bolus = strjoin(arrayfun(@char,this.bolus(:,ievent),'UniformOutput',false),',');
        fprintf(fid,['model.event(' num2str(ievent) ') = amievent(' ...
            str_trigger ', ...\n' ...
            '[' str_bolus '], ...\n' ...
            '[]);\n']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    fprintf(fid,'function r = pow(x,y)\n');
    fprintf(fid,'\n');
    fprintf(fid,'    r = x^y;\n');
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    fprintf(fid,'function r = power(x,y)\n');
    fprintf(fid,'\n');
    fprintf(fid,'    r = x^y;\n');
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    
    for ifun = 1:length(this.funmath)
        fprintf(fid,['function r = ' this.funarg{ifun} '\n']);
        fprintf(fid,'\n');
        fprintf(fid,['    r = ' this.funmath{ifun} ';\n']);
        fprintf(fid,'\n');
        fprintf(fid,'end\n');
        fprintf(fid,'\n');
    end
    
    fprintf(fid,['function r = factorial(x)\n']);
    fprintf(fid,'\n');
    fprintf(fid,['    error(''AMICI does not support the factorial function'');\n']);
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    
    fprintf(fid,['function r = ceil(x)\n']);
    fprintf(fid,'\n');
    fprintf(fid,['    error(''AMICI does not support the ceil function'');\n']);
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    
    fprintf(fid,['function r = psi(x)\n']);
    fprintf(fid,'\n');
    fprintf(fid,['    error(''AMICI does not support the psi function'');\n']);
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    
    fclose(fid);
end

function writeDefinition(header,identifier,field,this,fid)
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['%%%%\n%% ' header '\n']);
    if(length(this.(field))>0)
        fprintf(fid,['syms ' strjoin(cellfun(@char,num2cell(this.(field)),'UniformOutput',false)) '\n']);
    end
    fprintf(fid,['model.sym.' identifier ' = [' strjoin(cellfun(@char,num2cell(this.(field)),'UniformOutput',false),',') '];\n']);
end

function writeDerived(header,identifier,field,this,fid)
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['%%%%\n%% ' header '\n']);
    fprintf(fid,'\n');
    if(strcmp(header,'OBSERVABLES'))
        fprintf(fid,['%% ' strjoin(cellfun(@char,num2cell(this.observable_name),'UniformOutput',false),'\n%% ')  '\n']);
    end
    fprintf(fid,['model.sym.' identifier ' = [' strjoin(cellfun(@char,num2cell(this.(field)),'UniformOutput',false),', ...\n') '];']);
end