function this = gccode(this,csym,funstr,cvar,fid)
% gccode transforms symbolic expressions into c code and writes the
% respective expression into a specified file
%
% Parameters:
%  funstr: function for which C code is to be generated @type string
%  csym: symbolic expression to be transform @type symbolic 
%  cvar: name of the assigned variable in C @type string
%  fid: file id in which the expression should be written @type fileid
%
% Return values:
%  this: model definition object @type amimodel


lcsym = length(csym);
ff = find(csym~=0);
for j=transpose(ff(:));
    % temp = sprintf(mupadmex('generate::C', [cvar '[' num2str(j) '] =
    % ' char(csym(j))], 0)); this is much too slow for large systems
    if( lcsym == 1)
        cstr = [cvar ' = ' char(vpa(csym)) ';'];
%         tmp_str = ccode(csym(j));
%         cstr = [ cvar tmp_str(5:end) ';'];
    else
        cstr = [cvar '[' num2str(j-1) '] = ' char(vpa(csym(j))) ';'];
%         tmp_str = ccode(csym(j));
%         cstr = [ cvar '[' num2str(j-1) '] ' tmp_str(5:end) ';'];
    end

    
    % parse pow(...,...)
    idx = strfind(cstr,'^');
    if(~isempty(idx))
        bro = regexp(cstr,'\(');%find bracket opening
        brc = regexp(cstr,'\)');%find bracket closing
        % compute nested bracket level of expression
        brl = zeros(length(cstr),1);
        l = 0;
        for k = 1:length(cstr)
            if(~isempty(find(k==bro,1)))
                l = l+1;
            end
            brl(k) = l;
            if(~isempty(find(k==brc,1)))
                l = l-1;
            end
        end
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
    cstr = regexprep(cstr,'D\(\[([0-9]*)\]\, ([\w]*)\)\(','D$2\($1,'); % fix derivatives
    cstr = regexprep(cstr,'dirac\(([1-2]*),','Ddirac\($1,'); % fix derivatives
    cstr = regexprep(cstr,'dirac\(([^,]*), ([1-2]*)\)','Ddirac\($2,$1\)'); % fix matlab <= R2014a
    cstr = regexprep(cstr,'abs\(','fabs\('); % fix abs->fabs

    if(~(length(cstr)==1 && isempty(cstr{1})))
        
        cstr = strrep(cstr, [cvar ' ='],[cvar '[0] =']);
        
        if(strcmp(cvar,'qBdot'))
            cstr = regexprep(cstr,'qBdot\[([0-9]*)\]','qBdot\[ip]');
        elseif(strcmp(cvar,'y') || strcmp(cvar,'dydp'))
            cstr = regexprep(cstr,'x\[([0-9]*)\]','x\[it+nt*$1\]');
            cstr = regexprep(cstr,'y\[([0-9]*)\]','y\[it+nt*$1\]');
            if(strcmp(cvar,'dydp'))
                cstr = regexprep(cstr,'dydp\[([0-9]*)\]','dydp\[$1+ip*ny\]');
            end
        elseif(strcmp(cvar,'sy'))
            cstr = regexprep(cstr,'sx\[([0-9]*)\]','sx\[it+nt*\($1+ip*nx\)\]');
            cstr = regexprep(cstr,'sy\[([0-9]*)\]','sy\[it+nt*\($1+ip*ny\)\]');
            % this needs to stay after the sx replacement, it cannot be put
            % together with the replacements in y and dydp
            cstr = regexprep(cstr,'x\[([0-9]*)\]','x\[it+nt*$1\]');
            cstr = regexprep(cstr,'y\[([0-9]*)\]','y\[it+nt*$1\]');
        elseif(strcmp(cvar,'dxdotdp'))
            cstr = regexprep(cstr, 'dxdotdp\[([0-9]*)\]', 'dxdotdp[\($1+ip*nx\)\]');
        elseif(strcmp(cvar,'sdeltadisc'))
            cstr = regexprep(cstr,'sdeltadisc\[([0-9]*)\]','sdeltadisc\[plist[ip]+np*$1\]');
        elseif(strcmp(cvar,'Jdata'))
            cstr = strrep(cstr, 'Jdata[', 'J->data[');
        elseif(strcmp(cvar,'JBdata'))
            cstr = strrep(cstr, 'JBdata[', 'JB->data[');
        elseif(strcmp(cvar,'Jv'))
            cstr = strrep(cstr, 'v[', 'v_tmp[');
        elseif(strcmp(cvar,'JvB'))
            cstr = strrep(cstr, 'v[', 'vB_tmp[');
        elseif(strcmp(cvar,'sroot') || strcmp(cvar,'srootval'))
            cstr = regexprep(cstr, [cvar '\[([0-9]*)\] = '], [ cvar '[nroots + nmaxroot*(ip*nr + $1)] = ']);
        elseif(strcmp(cvar,'dsigma_tdp')  || strcmp(cvar,'dtdp')|| strcmp(cvar,'drvaldp'))
            cstr = regexprep(cstr, [cvar '\[([0-9]*)\] = '], [ cvar '[ip*nr + $1] = ']);
        elseif(strcmp(cvar,'dsigma_ydp'))
            cstr = regexprep(cstr, [cvar '\[([0-9]*)\] = '], [ cvar '[ip*ny + $1] = ']);
        end

        if(strfind(this.getFunArgs(funstr),'N_Vector'))
            cstr = strrep(cstr, 'dot[', 'dot_tmp[');
            cstr = strrep(cstr, 'x[', 'x_tmp[');
            cstr = strrep(cstr, 'B[', 'B_tmp[');
            cstr = strrep(cstr, '0[', '0_tmp[');
        end
        
    end
    fprintf(fid,[cstr '\n']);
end
end