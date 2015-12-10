function writeCcode(this,model, fid)
% writeCcode is a wrapper for gccode which initialises data and reduces
% overhead by check nonzero values
%
% Parameters: 
%  model: model defintion object @type amimodel
%  fid: file id in which the final expression is written @type fileid
%
%

ndisc = model.ndisc;
if(strcmp(funstr,'JSparse'))
    tmpfun = this;
    tmpfun.sym = this.sym(model.sparseidx);
    tmpfun.gccode(fid);
elseif(strcmp(funstr,'JSparseB'))
    tmpfun = this;
    tmpfun.sym = this.sym(model.sparseidxB);
    tmpfun.gccode(fid);
elseif(strcmp(funstr,'deltadisc') || strcmp(funstr,'sdeltadisc') || strcmp(funstr,'bdeltadisc') || strcmp(funstr,'ideltadisc'))
    nonzero = this.sym ~=0;
    if(any(any(nonzero)))
        fprintf(fid,'              switch(idisc) { \n');
        tmpfun = this;
        for idisc=1:ndisc
            if(any(nonzero(:,idisc)~=0))
                fprintf(fid,['              case ' num2str(idisc-1) ': {\n']);
                tmpfun.sym = this.sym(:,idisc);
                tmpfun.gccode(fid);
                fprintf(fid,'\n');
                fprintf(fid,'              } break;\n\n');
            end
        end
        fprintf(fid,'              } \n');
    end
else
    this.gccode(fid)
end


end