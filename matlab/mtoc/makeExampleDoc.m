opts.imageFormat = 'png';

examples = {'adjoint','dirac','dirac_adjoint','dirac_secondorder','dirac_secondorder_vectmult','events','jakstat_adjoint','steadystate'};

for ie = 1:length(examples)
    wrap_path = fileparts(which('amiwrap.m'));
    opts.format = 'html';
    publish(fullfile(wrap_path,'examples',['example_' examples{ie}],['example_' examples{ie} '.m']),opts)
    opts.format = 'html';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' examples{ie}],['model_' examples{ie} '_syms.m']),opts)
    opts.format = 'latex';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' examples{ie}],['example_' examples{ie} '.m']),opts)
    opts.format = 'latex';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' examples{ie}],['model_' examples{ie} '_syms.m']),opts)
end

wrap_path = fileparts(which('amiwrap.m'));

for ie = 1:length(examples)
    clear A
    fid = fopen(fullfile(wrap_path,'examples',['example_' examples{ie}],'html',['example_' examples{ie} '.html']),'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        str = tline;
        str = strrep(str,'src="example',['src="../examples/example_' examples{ie} '/html/example']);
        A{i} = str;
    end
    fclose(fid);
    
    % throw out file index
    
    % write changed file
    fid = fopen(fullfile(wrap_path,'examples',['example_' examples{ie}],'html',['example_' examples{ie} '.html']), 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);
end

for ie = 1:length(examples)
    clear A
    fid = fopen(fullfile(wrap_path,'examples',['example_' examples{ie}],'html',['example_' examples{ie} '.tex']),'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        str = tline;
        str = strrep(str,'\color{lightgray}','');
        str = strrep(str,'\color{black}','');
        str = strrep(str,'\begin{verbatim}','\begin{DoxyCode}');
        str = strrep(str,'\end{verbatim}','\end{DoxyCode}');
        str = strrep(str,'\includegraphics [width=4in]{example',['\includegraphics[width=\textwidth]{../../examples/example_' examples{ie} '/html/example']);
        if(all([~strfind(tline,'documentclass'),~strfind(tline,'usepackage'),~strfind(tline,'{document}'),~strfind(tline,'\sloppy'),~strfind(tline,'\setlength'),~strfind(tline,'\definecolor')]))
            A{i} = str;
        end
    end
    fclose(fid);
    
    % throw out file index
    
    % write changed file
    fid = fopen(fullfile(wrap_path,'examples',['example_' examples{ie}],'html',['example_' examples{ie} '.tex']), 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);
end

for ie = 1:length(examples)
    
    clear A
    
    fid = fopen(fullfile(wrap_path,'examples',['example_' examples{ie}],'html',['model_' examples{ie} '_syms.tex']),'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        str = tline;
        str = strrep(str,'\color{lightgray}','');
        str = strrep(str,'\color{black}','');
        str = strrep(str,'\begin{verbatim}','\begin{DoxyCode}');
        str = strrep(str,'\end{verbatim}','\end{DoxyCode}');
        if(all([~strfind(tline,'documentclass'),~strfind(tline,'usepackage'),~strfind(tline,'{document}'),~strfind(tline,'\sloppy'),~strfind(tline,'\setlength'),~strfind(tline,'\definecolor')]))
            A{i} = str;
        end
    end
    fclose(fid);
    
    % throw out file index
    
    % write changed file
    fid = fopen(fullfile(wrap_path,'examples',['example_' examples{ie}],'html',['model_' examples{ie} '_syms.tex']), 'w');
    for i = 1:numel(A)
        if A{i+1} == -1
            fprintf(fid,'%s', A{i});
            break
        else
            fprintf(fid,'%s\n', A{i});
        end
    end
    fclose(fid);
end

