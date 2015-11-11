opts.imageFormat = 'png';
for ie = 1:1
    save('ie.mat','ie','opts')
    wrap_path = fileparts(which('amiwrap.m'));
    opts.format = 'html';
    publish(fullfile(wrap_path,'examples',['example_' num2str(ie)],['example_model_' num2str(ie) '.m']),opts)
    load('ie.mat')
    opts.format = 'html';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' num2str(ie)],['example_model_' num2str(ie) '_syms.m']),opts)
    load('ie.mat')
    opts.format = 'latex';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' num2str(ie)],['example_model_' num2str(ie) '.m']),opts)
    load('ie.mat')
    opts.format = 'latex';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' num2str(ie)],['example_model_' num2str(ie) '_syms.m']),opts)
    load('ie.mat')
end

wrap_path = fileparts(which('amiwrap.m'));

for ie = 1:6
    clear A
    fid = fopen(fullfile(wrap_path,'examples',['example_' num2str(ie)],'html',['example_model_' num2str(ie) '.html']),'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        str = tline;
        str = strrep(str,'src="example',['src="../examples/example_' num2str(ie) '/html/example']);
        A{i} = str;
    end
    fclose(fid);
    
    % throw out file index
    
    % write changed file
    fid = fopen(fullfile(wrap_path,'examples',['example_' num2str(ie)],'html',['example_model_' num2str(ie) '.html']), 'w');
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

for ie = 1:6
    clear A
    fid = fopen(fullfile(wrap_path,'examples',['example_' num2str(ie)],'html',['example_model_' num2str(ie) '.tex']),'r');
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
        str = strrep(str,'\includegraphics [width=4in]{example',['\includegraphics[width=\textwidth]{../../examples/example_' num2str(ie) '/html/example']);
        if(all([~strfind(tline,'documentclass'),~strfind(tline,'usepackage'),~strfind(tline,'{document}'),~strfind(tline,'\sloppy'),~strfind(tline,'\setlength'),~strfind(tline,'\definecolor')]))
            A{i} = str;
        end
    end
    fclose(fid);
    
    % throw out file index
    
    % write changed file
    fid = fopen(fullfile(wrap_path,'examples',['example_' num2str(ie)],'html',['example_model_' num2str(ie) '.tex']), 'w');
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

for ie = 1:6
    
    clear A
    
    fid = fopen(fullfile(wrap_path,'examples',['example_' num2str(ie)],'html',['example_model_' num2str(ie) '_syms.tex']),'r');
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
    fid = fopen(fullfile(wrap_path,'examples',['example_' num2str(ie)],'html',['example_model_' num2str(ie) '_syms.tex']), 'w');
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

