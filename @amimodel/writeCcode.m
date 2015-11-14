function this = writeCcode(this, funstr, fid, ip, jp)
% writeCcode is a wrapper for this.gccode which initialises data and reduces
% overhead by check nonzero values
%
% Parameters: 
%  funstr: function for which C code is to be generated @type string
%  fid: file id in which the final expression is written @type fileid
%  ip: index for symbolic variable, if applicable svar(:,ip) will be passed to ggcode
%  jp: index for symbolic variable, if applicable svar(:,ip,jp) will be passed to ggcode
%
% Return values:
%  this: model definition object @type amimodel
%

ndisc = this.ndisc;
if(strcmp(funstr,'JSparse'))
    cvar =  'Jdata';
    this.gccode(this.sym.J(this.sparseidx),funstr,cvar,fid);
elseif(strcmp(funstr,'JSparseB'))
    cvar =  'JBdata';
    this.gccode(this.sym.JB(this.sparseidxB),funstr,cvar,fid);
elseif(strcmp(funstr,'deltadisc'))
    cvar =  'deltadisc';
    nonzero = this.sym.deltadisc ~=0;
    if(any(any(nonzero)))
        fprintf(fid,'              switch(idisc) { \n');
        for idisc=1:ndisc
            if(any(nonzero(:,idisc)~=0))
                fprintf(fid,['              case ' num2str(idisc-1) ': {\n']);
                this.gccode(this.sym.deltadisc(:,idisc),funstr,cvar,fid);
                fprintf(fid,'\n');
                fprintf(fid,'              } break;\n\n');
            end
        end
        fprintf(fid,'              } \n');
    end
elseif(strcmp(funstr,'sdeltadisc'))
    cvar =  'sdeltadisc';
    nonzero = this.sym.sdeltadisc(:,ip,:) ~=0;
    if(any(any(nonzero)))
        fprintf(fid,'              switch(idisc) { \n');
        for idisc=1:ndisc
            if(any(nonzero(:,idisc)~=0))
                fprintf(fid,['              case ' num2str(idisc-1) ': {\n']);
                this.gccode(this.sym.sdeltadisc(:,ip,idisc),funstr,cvar,fid);
                fprintf(fid,'\n');
                fprintf(fid,'              } break;\n\n');
            end
        end
        fprintf(fid,'              } \n');
    end
elseif(strcmp(funstr,'bdeltadisc'))
    cvar =  'bdeltadisc';
    nonzero = this.sym.bdeltadisc ~=0;
    if(any(any(nonzero)))
        fprintf(fid,'              switch(idisc) { \n');
        for idisc=1:ndisc
            if(any(nonzero(:,idisc)~=0))
                fprintf(fid,['              case ' num2str(idisc-1) ': {\n']);
                this.gccode(this.sym.bdeltadisc(:,idisc),funstr,cvar,fid);
                fprintf(fid,'\n');
                fprintf(fid,'              } break;\n\n');
            end
        end
        fprintf(fid,'              } \n');
    end
elseif(strcmp(funstr,'ideltadisc'))
    cvar =  'ideltadisc';
    nonzero = this.sym.ideltadisc ~=0;
    if(any(any(nonzero)))
        fprintf(fid,'              switch(idisc) { \n');
        for idisc=1:ndisc
            if(any(nonzero(:,idisc)~=0))
                fprintf(fid,['              case ' num2str(idisc-1) ': {\n']);
                this.gccode(this.sym.ideltadisc(:,idisc),funstr,cvar,fid);
                fprintf(fid,'\n');
                fprintf(fid,'              } break;\n\n');
            end
        end
        fprintf(fid,'              } \n');
    end
else
    if(nargin > 3)
        S.type = '()';
        if(nargin > 4)
            S.subs = {':',ip,jp};
        else
             S.subs = {':',ip};
        end
        cvar = this.getFunCVar(funstr);
        svar = this.getFunSVar(funstr);
        this.gccode(subsref(this.sym.(svar),S),funstr,cvar,fid);
    else
        cvar = this.getFunCVar(funstr);
        svar = this.getFunSVar(funstr);
        if(any(any(this.sym.(svar)~=0)))
            this.gccode(this.sym.(svar),funstr,cvar,fid);
        end
    end
end


end