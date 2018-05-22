function writeCcode_sensi(this,model,fid)
% writeCcode_sensi is a wrapper for writeCcode which loops over parameters and reduces
% overhead by check nonzero values
%
% Parameters:
%  model: model defintion object @type amimodel
%  fid: file id in which the final expression is written @type fileid
% 
% Return values:
%  void
      
np = model.np;
ng = model.ng;

nonzero_idx = find(this.sym);
nonzero = zeros(size(this.sym));
nonzero(nonzero_idx) = 1;

if(strcmp(this.funstr,'deltaqB'))
    if(any(nonzero))
        tmpfun = this;
        for ip=1:np*ng
            if(nonzero(ip))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                tmpfun.sym = this.sym(ip,:);
                tmpfun.writeCcode(model,fid);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(this.funstr,'deltasx'))
    if(any(any(any(nonzero))))
        tmpfun = this;
        for ip=1:np
            if(any(any(any(nonzero(:,ip,:)))))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                tmpfun.sym = permute(this.sym(:,ip,:),[1,3,2]);
                tmpfun.writeCcode(model,fid);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(this.funstr,'s2root'))
    if(any(any(any(nonzero))))
        tmpfun = this;
        for ip=1:np
            if(any(any(any(nonzero(:,ip,:)))))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                tmpfun.sym = squeeze(this.sym(model.z2event(1:model.nztrue),ip,:));
                tmpfun.writeCcode(model,fid);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(this.funstr,'qBdot'))
    nonzero = this.sym ~=0;
    if(any(any(nonzero)))
        tmpfun = this;
        for ip=1:np
            if(any(nonzero(:,ip)))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                tmpfun.sym = this.sym(:,ip);
                tmpfun.writeCcode(model,fid);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
else 
    nonzero = this.sym ~=0;
    if(any(any(nonzero)))
        tmpfun = this;
        for ip=1:np
            if(any(nonzero(:,ip)))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                if(strcmp(this.funstr,'sroot'))
                    tmpfun.sym = this.sym(model.z2event,ip);
                else
                    tmpfun.sym = this.sym(:,ip);
                end
                tmpfun.writeCcode(model,fid);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
end
end