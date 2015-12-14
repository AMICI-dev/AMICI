function writeCcode(this,model, fid)
% writeCcode is a wrapper for gccode which initialises data and reduces
% overhead by check nonzero values
%
% Parameters: 
%  model: model defintion object @type amimodel
%  fid: file id in which the final expression is written @type fileid
%
%

nevent = model.nevent;
if(strcmp(this.funstr,'JSparse'))
    tmpfun = this;
    tmpfun.sym = model.fun.J.sym(model.sparseidx);
    tmpfun.gccode(fid);
elseif(strcmp(this.funstr,'JSparseB'))
    tmpfun = this;
    tmpfun.sym = model.fun.JB.sym(model.sparseidxB);
    tmpfun.gccode(fid);
elseif(strcmp(this.funstr,'deltadisc') || strcmp(this.funstr,'sdeltadisc') || strcmp(this.funstr,'bdeltadisc') || strcmp(this.funstr,'ideltadisc'))
    nonzero = this.sym ~=0;
    if(any(any(nonzero)))
        fprintf(fid,'              switch(ie) { \n');
        tmpfun = this;
        for ie=1:nevent
            if(any(nonzero(:,ie)~=0))
                fprintf(fid,['              case ' num2str(ie-1) ': {\n']);
                tmpfun.sym = this.sym(:,ie);
                tmpfun.gccode(fid);
                fprintf(fid,'\n');
                fprintf(fid,'              } break;\n\n');
            end
        end
        fprintf(fid,'              } \n');
    end
else
    this.gccode(fid);
end


end