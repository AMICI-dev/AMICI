function this = gccode(this,fid)
% gccode transforms symbolic expressions into c code and writes the
% respective expression into a specified file
%
% Parameters:
%  fid: file id in which the expression should be written @type fileid
%
% Return values:
%  this: function definition object @type amifun

ff = find(this.sym~=0);
for j=transpose(ff(:));
    % temp = sprintf(mupadmex('generate::C', [cvar '[' num2str(j) '] =
    % ' char(csym(j))], 0)); this is much too slow for large systems so we
    % built our own implementation here
    
    %% 
    % convert to strings
    
    % apply vpa to convert to fix point numbers, this prevents
    % erroneous computations from bad casting.
    % we always output vectors!
    cstr = [this.cvar '[' num2str(j-1) '] = ' char(vpa(this.sym(j))) ';'];
    
    %%
    % parse pow(...,...)
    
    % find ^
    idx = strfind(cstr,'^');
    % parse brackets
    if(~isempty(idx))
        % compute bracket levels add one for each (, remove 1 for each )
        brl = cumsum(cstr == '(') - cumsum(cstr == ')');
    end
    while(~isempty(idx))
        samebrl = find(brl == brl(idx(1))); % find bracket level same as found level
        % find first argument
        if(strcmp(cstr(idx(1)-1),')'))% if the expression is braced find start of higher levels
            arg1start = samebrl(find(samebrl<idx(1),1,'last')) + 1;
            while(~isempty(regexp(cstr(arg1start-1),'[\w]','once')))
                arg1start = arg1start-1;
            end
        else % otherwise find start of expression
            m = regexp(cstr(1:(idx(1)-1)),'([\w\[\]]*)','start');
            arg1start = m(end);
        end
        arg1end = idx(1)-1;
        % find second argument
        arg2start = idx(1)+1;
        if(strcmp(cstr(idx(1)+1),'('))% if the expression is braced find end of higher levels
            arg2end = samebrl(find(samebrl>idx(1),1,'first')) - 1;
        else % otherwise find end of expression
            m = regexp(cstr((idx(1)+1):end),'([\w\[\]]*)','end');
            arg2end = idx(1)+m(1);
        end
        
        % recombine
        cstr = [cstr(1:(arg1start-1)) '(pow('                                                    cstr(arg1start:arg1end)   ','              cstr(arg2start:arg2end)                 '))'             cstr((arg2end+1):end)];
        % compute new bracket level
        brl  = [ brl(1:(arg1start-1)); brl(idx(1))+1; ones(3,1)*brl(idx(1))+1 ; brl(idx(1))+2 ; brl(arg1start:arg1end)+2 ; brl(idx(1))+2 ; brl(arg2start:arg2end)+2 ; brl(idx(1))+2; brl(idx(1))+1; brl((arg2end+1):end) ];
        idx = strfind(cstr,'^');
    end
    cstr = deblank(cstr);
%     cstr = regexprep(cstr,'D\(\[([0-9]*)\]\, ([\w]*)\)\(','D$2\($1,'); % fix derivatives
%     cstr = regexprep(cstr,'dirac\(([1-2]*),','Ddirac\($1,'); % fix derivatives
%     cstr = regexprep(cstr,'dirac\(([^,]*), ([1-2]*)\)','Ddirac\($2,$1\)'); % fix matlab <= R2014a
    cstr = regexprep(cstr,'abs\(','fabs\('); % fix abs->fabs

    if(~(length(cstr)==1 && isempty(cstr{1})))
        
        % fix various function specific variable names/indexes
        
        cstr = regexprep(cstr,'x_([0-9]+)','x\[$1\]');
        cstr = regexprep(cstr,'p_([0-9]+)','p\[$1\]');
        cstr = regexprep(cstr,'k_([0-9]+)','k\[$1\]');
        cstr = regexprep(cstr,'xB_([0-9]+)','xB\[$1\]');
        cstr = regexprep(cstr,'J([0-9]+)','J\[$1\]');
        cstr = regexprep(cstr,'xdotdp([0-9]+)','xdotdp\[$1 + ip*nx\]');
        
        if(strcmp(this.cvar,'qBdot'))
            cstr = regexprep(cstr,'qBdot\[([0-9]*)\]','qBdot\[ip]');
        elseif(strcmp(this.cvar,'y') || strcmp(this.cvar,'dydp'))
            cstr = regexprep(cstr,[this.cvar '\[([0-9]*)\]'],[this.cvar '\[it+nt*$1\]']);
        elseif(strcmp(this.cvar,'sy') || strcmp(this.cvar,'dydp'))
            cstr = regexprep(cstr,[this.cvar '\[([0-9]*)\]'],[this.cvar '\[it+nt*\($1+ip*ny\)\]']);
        elseif(strcmp(this.cvar,'dxdotdp'))
            cstr = regexprep(cstr, 'dxdotdp\[([0-9]*)\]', 'dxdotdp[\($1+ip*nx\)\]');
        elseif(strcmp(this.cvar,'sdeltadisc'))
            cstr = regexprep(cstr,'sdeltadisc\[([0-9]*)\]','sdeltadisc\[plist[ip]+np*$1\]');
        elseif(strcmp(this.cvar,'Jdata'))
            cstr = strrep(cstr, 'Jdata[', 'J->data[');
        elseif(strcmp(this.cvar,'JBdata'))
            cstr = strrep(cstr, 'JBdata[', 'JB->data[');
        elseif(strcmp(this.cvar,'Jv'))
            cstr = strrep(cstr, 'v[', 'v_tmp[');
        elseif(strcmp(this.cvar,'JvB'))
            cstr = strrep(cstr, 'v[', 'vB_tmp[');
        elseif(strcmp(this.cvar,'sxdot'))
            cstr = strrep(cstr,'tmp_J[','tmp_J->data[');
        elseif(strcmp(this.cvar,'sroot') || strcmp(this.cvar,'srootval'))
            cstr = regexprep(cstr, [this.cvar '\[([0-9]*)\] = '], [ this.cvar '[nroots + nmaxroot*(ip*nr + $1)] = ']);
        elseif(strcmp(this.cvar,'dsigma_zdp')  || strcmp(this.cvar,'sz')|| strcmp(this.cvar,'dzdp'))
            cstr = regexprep(cstr, [this.cvar '\[([0-9]*)\] = '], [ this.cvar '[ip*ne + $1] = ']);
        elseif(strcmp(this.cvar,'dsigma_ydp'))
            cstr = regexprep(cstr, [this.cvar '\[([0-9]*)\] = '], [ this.cvar '[ip*ny + $1] = ']);
        end

        
        for nvec = this.nvecs
            cstr = strrep(cstr, [nvec{1} '['],[nvec{1} '_tmp[']);
        end
        if(strcmp(this.funstr,'dydx'))
            cstr = strrep(cstr,'dydx_tmp','dydx');
        end
        
    end
    % print to file
    fprintf(fid,[cstr '\n']);
end
end