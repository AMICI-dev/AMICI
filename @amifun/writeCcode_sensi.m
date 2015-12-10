function this = writeCcode_sensi(this,model,fid)
% writeCcode_sensi is a wrapper for writeCcode which loops over parameters and reduces
% overhead by check nonzero values
%
% Parameters:
%  svar: symbolic variable which is later on passed to ggode @type symbolic
%  fid: file id in which the final expression is written @type fileid
% 
% Return values:
%  this: model definition object @type amimodel
      
np = model.np;

nonzero = this.sym~=0;

if(strcmp(svar,'drvaldp'))
    
    if(any(any(nonzero)))
        tmpfun = this;
        for ip=1:np
            if(any(nonzero(:,ip)))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                tmpfun.sym = this.sym(:,ip);
                tmpfun.writeCcode(fid);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(svar,'s2root'))
    if(any(any(any(nonzero))))
        tmpfun = this;
        for ip=1:np
            if(any(any(nonzero(:,:,ip))))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                fprintf(fid,'    sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                fprintf(fid,['      switch(jp) {\n']);
                for jp=1:np
                    if(any(nonzero(:,jp,ip)))
                        fprintf(fid,['          case ' num2str(jp-1) ': {\n']);
                        tmpfun.sym = this.sym(:,ip,jp);
                        tmpfun.writeCcode(fid);
                        fprintf(fid,'\n');
                        fprintf(fid,'          } break;\n\n');
                    end
                end
                fprintf(fid,'\n      }\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(svar,'s2rootval'))
    if(any(any(any(nonzero))))
        for ip=1:np
            if(any(any(nonzero(:,:,ip))))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                fprintf(fid,'  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                fprintf(fid,['      switch(jp) {\n']);
                for jp=1:np
                    if(any(nonzero(:,jp,ip)))
                        fprintf(fid,['          case ' num2str(jp-1) ': {\n']);
                        this.writeCcode('s2rootval', fid, ip, jp);
                        fprintf(fid,'\n');
                        fprintf(fid,'           } break;\n\n');
                    end
                end
                fprintf(fid,'\n      }\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(svar,'ideltadisc'))
    if(any(nonzero))
        for ip=1:np
            if(nonzero(ip))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                this.writeCcode(svar, fid, ip);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
elseif(strcmp(svar,'sdeltadisc'))
    if(any(any(any(nonzero))))
        for ip=1:np
            if(any(any(any(nonzero(:,ip,:)))))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                if(strcmp(svar,'sroot') || strcmp(svar,'srootval'))
                    fprintf(fid,'  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                end
                this.writeCcode(svar, fid, ip);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
else 
    nonzero = this.sym.(svar) ~=0;
    if(any(any(nonzero)))
        for ip=1:np
            if(any(nonzero(:,ip)))
                fprintf(fid,['  case ' num2str(ip-1) ': {\n']);
                if(strcmp(svar,'sroot') || strcmp(svar,'srootval'))
                    fprintf(fid,'  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);\n');
                end
                this.writeCcode(svar, fid, ip);
                fprintf(fid,'\n');
                fprintf(fid,'  } break;\n\n');
            end
        end
    end
end
end