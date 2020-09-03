function writeCcode(this,model,fid)
% writeCcode is a wrapper for gccode which initialises data and reduces
% overhead by check nonzero values
%
% Parameters:
%  model: model defintion object @type amimodel
%  fid: file id in which the final expression is written @type fileid
%
% Return values:
%  void

nonzero_idx = find(this.sym);
nonzero = zeros(size(this.sym));
nonzero(nonzero_idx) = 1;

nevent = model.nevent;
if(strcmp(this.funstr,'JSparse'))
    tmpfun = this;
    tmpfun.sym = model.fun.J.sym(model.sparseidx);
    tmpfun.gccode(model,fid);
elseif(ismember(this.funstr,{'rz','z','sz','srz'}))
    if(any(nonzero))
        fprintf(fid,'    switch(ie) { \n');
        for ievent=1:nevent
            tmpfun = this;
            % set all z that do not belong to this event to zero
            % dont shorten the vector as we need the indices
            fprintf(fid,['        case ' num2str(ievent-1) ': {\n']);
            tmpfun.sym(model.z2event~=ievent) = 0;
            tmpfun.gccode(model,fid);
            fprintf(fid,'\n');
            fprintf(fid,'        } break;\n\n');
        end
        fprintf(fid,'    } \n');
    end
elseif(strcmp(this.funstr,'stau') )
    if(any(nonzero))
        fprintf(fid,'    switch(ie) { \n');
        for ievent=1:nevent
            tmpfun = this;
            % set all z that do not belong to this event to zero
            % dont shorten the vector as we need the indices
            fprintf(fid,['        case ' num2str(ievent-1) ': {\n']);
            tmpfun.sym = tmpfun.sym(ievent);
            tmpfun.gccode(model,fid);
            fprintf(fid,'\n');
            fprintf(fid,'        } break;\n\n');
        end
        fprintf(fid,'    } \n');
    end
elseif(strcmp(this.funstr,'deltax') || strcmp(this.funstr,'deltasx') || strcmp(this.funstr,'deltaxB') || strcmp(this.funstr,'deltaqB'))
    if(any(any(nonzero)))
        fprintf(fid,'              switch(ie) { \n');
        tmpfun = this;
        for ievent=1:nevent
            if(any(nonzero(:,ievent)))
                fprintf(fid,['              case ' num2str(ievent-1) ': {\n']);
                tmpfun.sym = this.sym(:,ievent);
                tmpfun.gccode(model,fid);
                fprintf(fid,'\n');
                fprintf(fid,'              } break;\n\n');
            end
        end
        fprintf(fid,'              } \n');
    end
elseif(any(strcmp(this.funstr,{'Jy','dJydsigma','dJydy'})))
    tmpfun = this;
    if(any(any(any(nonzero))))
        fprintf(fid,['switch(iy){\n']);
        for iy = 1:model.nytrue
            fprintf(fid,['    case ' num2str(iy-1) ':\n']);
            tmpfun.sym = permute(this.sym(iy,:,:),[2,3,1]);
            tmpfun.gccode(model,fid);
            fprintf(fid,'    break;\n');
        end
        fprintf(fid,['}\n']);
    end
elseif(any(strcmp(this.funstr,{'Jz','dJzdsigma','dJzdz','Jrz','dJrzdsigma','dJrzdz'})))
    tmpfun = this;
    if(any(any(any(nonzero))))
        fprintf(fid,['switch(iz){\n']);
        for iz = 1:model.nztrue
            fprintf(fid,['    case ' num2str(iz-1) ':\n']);
            tmpfun.sym = permute(this.sym(iz,:,:),[2,3,1]);
            tmpfun.gccode(model,fid);
            fprintf(fid,'    break;\n');
        end
        fprintf(fid,['}\n']);
    end
else
    if(any(any(any(nonzero))))
        this.gccode(model,fid);
    end
end


end