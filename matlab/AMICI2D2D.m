function AMICI2D2D( filename, modelname )
    
    eval(['model = ' filename '();'])
    
    fid = fopen([modelname '.def'],'w');
    
    fprintf(fid,'DESCRIPTION\n');
    fprintf(fid,['"' modelname '"\n']);
    fprintf(fid,'\n');
    fprintf(fid,'PREDICTOR\n');
    fprintf(fid,'t               T   min         time	0	60');
    fprintf(fid,'\n');
    fprintf(fid,'COMPARTMENTS\n');
    fprintf(fid,'cell             V   "pl"      "vol."        1\n');
    fprintf(fid,'\n');
    fprintf(fid,'STATES\n');
    fprintf(fid,strjoin(cellfun(@(x) [char(x) '\t\tC\t"nM"\t"conc."\tcell\n'],num2cell(model.sym.x),'UniformOutput',false),''));
    fprintf(fid,'\n');
    fprintf(fid,'INPUTS\n');
    fprintf(fid,'\n');
    fprintf(fid,'ODES\n');
    tmp = strjoin(cellfun(@(x) ['"' char(x) '"\n'],num2cell(model.sym.xdot),'UniformOutput',false),'');
    tmp = regexprep(tmp,'heaviside\(t - ([0-9/]*)\)','step1(t,0.0,$1,1.0)');
    tmp = regexprep(tmp,'heaviside\(([0-9/]*) - t\)','step1(t,1.0,$1,0.0)');
    fprintf(fid,tmp);
    fprintf(fid,'\n');
    fprintf(fid,'DERIVED\n');
    fprintf(fid,'\n');
    fprintf(fid,'OBSERVABLES\n');
    fprintf(fid,strjoin(cellfun(@(x,y) ['y_' num2str(y) '\tC\t"au"\t"conc."\t1\t0\t"' char(x) '"\n'],num2cell(model.sym.y),num2cell(1:length(model.sym.y)),'UniformOutput',false),''));
    fprintf(fid,'\n');
    fprintf(fid,'ERRORS\n');
    fprintf(fid,strjoin(cellfun(@(x) ['y_' num2str(x) '\t"sd_y_' num2str(x) '"\n'],num2cell(1:length(model.sym.y)),'UniformOutput',false),''));
    fprintf(fid,'\n');
    fprintf(fid,'CONDITIONS\n');
    fprintf(fid,strjoin(cellfun(@(x,y) ['init_' char(x) '\t"' char(sym(y)) '"\n'],num2cell(model.sym.x),num2cell(model.sym.x0),'UniformOutput',false),''));
    fprintf(fid,strjoin(cellfun(@(x) [ char(x) '\t"1"\n'],num2cell(model.sym.k),'UniformOutput',false),''));
    fclose(fid);
end