function writeAMICI(this,modelname)
    % importSBML writes the information stored in the SBMLode object into
    % an AMICI model definition file
    %
    % Parameters:
    %  modelname: target name of the AMICI model definition file, _syms
    %  will be appended
    %
    % Return values:
    %  this: initialised SBMLode object object @type SBMLode
    fprintf('writing file ...\n')
    fid = fopen([modelname '_syms.m'],'w');
    
    fprintf(fid,['function model = ' modelname '_syms()\n']);
    fprintf(fid,'\n');
    if(strcmp(this.time_symbol,''))
        fprintf(fid,'t = sym(''t'');\n');
    else
        fprintf(fid,[this.time_symbol ' = sym(''t'');\n']);
    end
    fprintf(fid,'\n');
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
        fprintf(fid,['model.sym.event(' num2str(ievent) ') = amivent(' ...
            char(this.trigger(ievent)) ', ...\n' ...
            '[' strjoin(cellfun(@char,{this.bolus(:,ievent)},'UniformOutput',false),',') '], ...\n' ...
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
    fclose(fid);
end

function writeDefinition(header,identifier,field,this,fid)
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['%%%%\n%% ' header '\n']);
    fprintf(fid,['syms ' strjoin(cellfun(@char,num2cell(this.(field)),'UniformOutput',false)) '\n']);
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