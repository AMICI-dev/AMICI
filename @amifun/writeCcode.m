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

nevent = model.nevent;
if(strcmp(this.funstr,'JSparse'))
    tmpfun = this;
    tmpfun.sym = model.fun.J.sym(model.sparseidx);
    tmpfun.gccode(model,fid);
elseif(strcmp(this.funstr,'JSparseB'))
    tmpfun = this;
    tmpfun.sym = model.fun.JB.sym(model.sparseidxB);
    tmpfun.gccode(model,fid);
elseif(strcmp(this.funstr,'z') || strcmp(this.funstr,'sz'))
    nonzero = this.sym ~=0;
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
elseif(strcmp(this.funstr,'deltax') || strcmp(this.funstr,'deltasx') || strcmp(this.funstr,'deltaxB') || strcmp(this.funstr,'deltaqB'))
    nonzero = this.sym ~=0;
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
elseif(strcmp(this.funstr,'Jy') || strcmp(this.funstr,'dJydp') || strcmp(this.funstr,'dJydx') || strcmp(this.funstr,'sJy'))
    nonzero = this.sym ~=0;
    if(any(any(nonzero)))
        fprintf(fid,'int iy;\n');
        fprintf(fid,'for(iy=0;iy<ny;iy++){\n');
        fprintf(fid,'    if(mxIsNaN(my[iy*nt+it])){\n');
        fprintf(fid,'        my[iy*nt+it] = y[iy*nt+it];\n');
        fprintf(fid,'    }\n');
        fprintf(fid,'}\n');
        this.gccode(fid);
    end
else
    nonzero = this.sym ~=0;
    if(any(any(nonzero)))
        this.gccode(model,fid);
    end
end


end